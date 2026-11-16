#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Ensembl Regulatory Build features overlapping a region (the skill's only own HTTP client).

Subcommands:
  overlap  GET /overlap/region/human/{region}?feature=regulatory and write the JSON list as served.

The ensembl-database skill has no Regulatory Build command, so this client follows that skill's
rule for custom queries: at most 15 requests per second (shared across processes through a file
lock), exponential backoff on 429/5xx and on network errors or timeouts (3 attempts), and the
response body in the error otherwise. Ensembl names the mitochondrion MT; chrM/M are mapped to it.
The output is written before the release lookup (/info/data), so a failure there only drops the
release from the success line.

Example:
  uv run scripts/regulatory_features.py overlap --region chr4:1000-2000 --output raw/regulatory.json
"""

import argparse
import json
import os
import re
import socket
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

try:
    import fcntl
except ImportError:  # pragma: no cover - non-POSIX: the lock degrades to per-process pacing
    fcntl = None

BASE_URL = "https://rest.ensembl.org"
LOCK_PATH = Path.home() / ".cache" / "alphagenome-gene-skills" / "ensembl.lock"
USER_AGENT = "alphagenome-gene-skills/regulatory_features (research use)"
_REGION_RE = re.compile(r"^(?:chr)?([0-9XYM]{1,2}|MT):(\d+)-(\d+)$", re.I)


class RateLimitError(RuntimeError):
    """Ensembl answered 429 Too Many Requests on every attempt."""


class EnsemblHttpError(RuntimeError):
    """Non-retryable HTTP error; the message carries the response body."""


class EnsemblRegulatoryClient:
    """Rate-limited Ensembl REST client for Regulatory Build overlap queries.

    Attributes:
      REQUESTS_PER_SECOND: Ensembl REST limit (ensembl-database skill: 15 requests/s).
      MAX_ATTEMPTS: Attempts per request on 429 or 5xx before giving up.
      TIMEOUT_SECONDS: Socket timeout per attempt.
    """

    REQUESTS_PER_SECOND = 15
    MAX_ATTEMPTS = 3
    TIMEOUT_SECONDS = 60
    BACKOFF_SECONDS = 1.0

    def __init__(self, base_url: str = BASE_URL, lock_path: Path = LOCK_PATH):
        self.base_url = base_url.rstrip("/")
        self.lock_path = lock_path

    def _pace(self) -> None:
        """Waits so that all processes together stay under REQUESTS_PER_SECOND."""
        interval = 1.0 / self.REQUESTS_PER_SECOND
        self.lock_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.lock_path, "a+", encoding="utf-8") as handle:
            if fcntl is not None:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                handle.seek(0)
                text = handle.read().strip()
                last = float(text) if text else 0.0
                wait = last + interval - time.time()
                if wait > 0:
                    time.sleep(wait)
                handle.seek(0)
                handle.truncate()
                handle.write(f"{time.time():.6f}")
                handle.flush()
            finally:
                if fcntl is not None:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    def get_json(self, path: str):
        """GETs `path` (relative to the base URL) and returns the decoded JSON."""
        url = f"{self.base_url}{path}"
        request = urllib.request.Request(
            url, headers={"Accept": "application/json", "User-Agent": USER_AGENT}
        )
        for attempt in range(1, self.MAX_ATTEMPTS + 1):
            self._pace()
            try:
                with urllib.request.urlopen(request, timeout=self.TIMEOUT_SECONDS) as response:
                    return json.loads(response.read().decode("utf-8"))
            except urllib.error.HTTPError as exc:
                body = exc.read().decode("utf-8", errors="replace")
                retryable = exc.code == 429 or 500 <= exc.code < 600
                if retryable and attempt < self.MAX_ATTEMPTS:
                    retry_after = exc.headers.get("Retry-After") if exc.headers else None
                    delay = float(retry_after) if retry_after and retry_after.replace(".", "", 1).isdigit() else 0.0
                    time.sleep(max(delay, self.BACKOFF_SECONDS * 2 ** (attempt - 1)))
                    continue
                if exc.code == 429:
                    raise RateLimitError(
                        f"HTTP 429 from Ensembl after {attempt} attempt(s) for {url}; wait and retry."
                    ) from exc
                raise EnsemblHttpError(f"HTTP {exc.code} from Ensembl for {url}: {body.strip()}") from exc
            except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError) as exc:
                if attempt < self.MAX_ATTEMPTS:
                    time.sleep(self.BACKOFF_SECONDS * 2 ** (attempt - 1))
                    continue
                reason = exc.reason if isinstance(exc, urllib.error.URLError) else exc
                raise EnsemblHttpError(f"Network error for {url} after {attempt} attempt(s): {reason}") from exc
        raise EnsemblHttpError(f"No response for {url}")  # pragma: no cover

    def release(self) -> int | None:
        """Ensembl release served by the REST API (/info/data), or None if not reported."""
        data = self.get_json("/info/data")
        releases = data.get("releases") if isinstance(data, dict) else None
        return int(max(releases)) if releases else None

    def overlap(self, region: str) -> list[dict]:
        """Regulatory Build features overlapping `region` (1-based closed, chr prefix optional)."""
        match = _REGION_RE.match(region.replace(",", "").strip())
        if not match:
            raise ValueError(f"Invalid region {region!r}; expected CHR:START-END (1-based closed).")
        chrom, start, end = match.group(1).upper(), int(match.group(2)), int(match.group(3))
        chrom = "MT" if chrom == "M" else chrom
        if end < start:
            raise ValueError(f"Invalid region {region!r}: end < start.")
        data = self.get_json(f"/overlap/region/human/{chrom}:{start}-{end}?feature=regulatory")
        if not isinstance(data, list):
            raise EnsemblHttpError(f"Unexpected response for {region}: expected a JSON list.")
        return data


def overlap(args: argparse.Namespace) -> None:
    client = EnsemblRegulatoryClient()
    features = client.overlap(args.region)
    output = Path(args.output)
    if output.parent and not output.parent.exists():
        os.makedirs(output.parent, exist_ok=True)
    output.write_text(json.dumps(features, indent=2) + "\n", encoding="utf-8")
    try:
        release = client.release()
    except (RateLimitError, EnsemblHttpError, ValueError) as exc:
        print(f"Warning: Ensembl release lookup failed ({exc}); record the release by hand.", file=sys.stderr)
        release = None
    print(f"Wrote {len(features)} regulatory feature(s) for {args.region} "
          f"(Ensembl release {release if release is not None else 'unknown'}) to {output}")


def main() -> None:
    for stream in (sys.stdout, sys.stderr):  # UTF-8 on any console (cp1252 on Windows)
        stream.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("overlap", help="Regulatory Build features overlapping a region")
    p.add_argument("--region", required=True, help="CHR:START-END, 1-based closed (e.g. chr4:1000-2000)")
    p.add_argument("--output", required=True, help="Output JSON path")
    p.set_defaults(func=overlap)

    args = parser.parse_args()
    try:
        args.func(args)
    except (RateLimitError, EnsemblHttpError, ValueError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

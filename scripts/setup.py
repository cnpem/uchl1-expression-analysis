#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Prepares a computer for the AlphaGenome analysis skills, with the researcher's approval for each item.

Standard library only, so it runs with `uv run` (uv brings its own Python) on Windows, macOS and Linux.
Nothing is installed without `--approved <item>` naming that same item: the agent asks the researcher
first, one item at a time, and passes the approval it received.

Subcommands:
  diagnose   What the computer has: system and architecture, uv, Node.js (npx), the analysis and
             dependency skills, terms of use accepted, keys (presence only), disk space, network.
  plan       What is missing, as questions (title, question, why, options with label and effect,
             recommended, handled_by), in the shape of `expression.py check --json`.
  terms      Shows the terms of use a dependency skill asks to be accepted (links), for the agent to
             show before asking.
  install    Installs or records one approved item: setup.uv, setup.node, setup.skill:<name>,
             setup.licence:<skill>, setup.key:<VAR> (creates ~/.env with an empty line and opens it in
             the system's text editor; the researcher pastes the key there; it never reaches the chat).
  key        Says whether a key is set (in the environment or ~/.env), never its value.
  verify     Runs each installed analysis skill's `check`, and its offline `reproduce` when the
             reference study is bundled.
  report     Writes environment-report.md in the workspace: what is installed, versions, the result
             of verify; for the researcher to keep or attach when asking for help.

Run from the researcher's workspace (the folder the agent works in). Writes only: skills installed by
`npx skills add` (into the agents' skill folders of the workspace), .licenses/ files the researcher
accepted, ~/.env (only when asked, with an empty key line), .alphagenome-setup/state.json and
environment-report.md. Never prints, logs or copies the value of a key.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import re
import shutil
import struct
import subprocess
import sys
import time
import urllib.error
import urllib.request
import webbrowser
from pathlib import Path
from typing import Any

SKILL_DIR = Path(__file__).resolve().parent.parent
TARGETS = json.loads((SKILL_DIR / "assets" / "targets.json").read_text(encoding="utf-8"))
STATE_DIR = ".alphagenome-setup"
REPORT_FILE = "environment-report.md"
UV_INSTALL = {
    "Windows": ["powershell", "-NoProfile", "-ExecutionPolicy", "ByPass", "-c", "irm https://astral.sh/uv/install.ps1 | iex"],
    "other": ["sh", "-c", "curl -LsSf https://astral.sh/uv/install.sh | sh"],
}
NODE_PAGE = "https://nodejs.org/en/download"
_LICENCE_RE = re.compile(r"\.licenses/([A-Za-z0-9_.-]+_LICENSE\.txt)")
_URL_RE = re.compile(r"https?://[^\s,)]+")


# ------------------------------------------------------------------------------------------- helpers


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _run(cmd: list[str], cwd: Path | None = None, timeout: int = 900) -> tuple[int, str]:
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, encoding="utf-8", errors="replace",
                           timeout=timeout)
        return p.returncode, (p.stdout or "") + (p.stderr or "")
    except (OSError, subprocess.TimeoutExpired) as exc:
        return 127, str(exc)


def _version(cmd: list[str]) -> str | None:
    code, out = _run(cmd, timeout=60)
    return out.strip().splitlines()[0] if code == 0 and out.strip() else None


def uv_path() -> str | None:
    found = shutil.which("uv")
    if found:
        return found
    for cand in (Path.home() / ".local" / "bin" / ("uv.exe" if os.name == "nt" else "uv"),
                 Path.home() / ".cargo" / "bin" / ("uv.exe" if os.name == "nt" else "uv")):
        if cand.is_file():
            return str(cand)
    return None


def skill_dirs(ws: Path) -> list[Path]:
    """Where agents look for skills: next to this skill, the workspace's agent folders, the user's."""
    home = Path.home()
    out = [SKILL_DIR.parent, ws / ".claude" / "skills", ws / ".agents" / "skills",
           home / ".claude" / "skills", home / ".agents" / "skills"]
    return list(dict.fromkeys(p.resolve() for p in out))


def find_skill(ws: Path, name: str) -> Path | None:
    for d in skill_dirs(ws):
        if (d / name / "SKILL.md").is_file():
            return d / name
    return None


def key_present(name: str) -> bool:
    """True when the key is set in the environment or in ~/.env. The value is never kept or shown."""
    if os.environ.get(name):
        return True
    try:
        for line in (Path.home() / ".env").read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip().removeprefix("export ").strip()
            key, sep, value = line.partition("=")
            if sep and key.strip() == name and value.strip().strip("\"'"):
                return True
    except OSError:
        pass
    return False


def terms_of(skill_dir: Path) -> tuple[str, list[str]] | None:
    """(licence file, links to the terms) a dependency skill asks for in its SKILL.md, or None."""
    text = (skill_dir / "SKILL.md").read_text(encoding="utf-8", errors="replace")
    m = _LICENCE_RE.search(text)
    if not m:
        return None
    para = text[m.start(): text.find("\n\n", m.start()) if "\n\n" in text[m.start():] else len(text)]
    links = [u.rstrip(".") for u in _URL_RE.findall(para)]
    return m.group(1), links


def arch() -> dict[str, str]:
    info = {"python_bits": f"{struct.calcsize('P') * 8}-bit", "python_machine": platform.machine()}
    if os.name == "nt":
        native = os.environ.get("PROCESSOR_ARCHITEW6432") or os.environ.get("PROCESSOR_ARCHITECTURE", "")
        ident = os.environ.get("PROCESSOR_IDENTIFIER", "")
        info["processor"] = ident
        info["native"] = "ARM64" if ("ARM" in ident.upper() or native.upper() == "ARM64") else native
    else:
        info["native"] = platform.machine()
    return info


def network(urls: list[str]) -> dict[str, str]:
    out = {}
    for u in urls:
        try:
            req = urllib.request.Request(u, headers={"User-Agent": "alphagenome-skills-setup"})
            with urllib.request.urlopen(req, timeout=10) as r:
                out[u] = f"ok ({r.status})"
        except urllib.error.HTTPError as exc:  # the server answered: reachable
            out[u] = f"ok ({exc.code})"
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            out[u] = f"unreachable ({getattr(exc, 'reason', exc)})"
    return out


def state(ws: Path) -> dict[str, Any]:
    p = ws / STATE_DIR / "state.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else {"installs": [], "verify": None}


def save_state(ws: Path, st: dict[str, Any]) -> None:
    (ws / STATE_DIR).mkdir(exist_ok=True)
    (ws / STATE_DIR / "state.json").write_text(json.dumps(st, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def targets_for(ws: Path, wanted: list[str] | None) -> list[str]:
    """Analysis skills to prepare: the ones named, else the ones installed, else expression-analysis."""
    names = list(TARGETS["analysis_skills"])
    if wanted:
        bad = [w for w in wanted if w not in names]
        if bad:
            sys.exit(f"Unknown analysis skill(s) {bad}; known: {names}")
        return wanted
    return [n for n in names if find_skill(ws, n)] or ["expression-analysis"]


# ------------------------------------------------------------------------------------------ diagnose


def diagnose_data(ws: Path, wanted: list[str] | None, net: bool = True) -> dict[str, Any]:
    targets = targets_for(ws, wanted)
    deps = list(dict.fromkeys(d for t in targets for d in TARGETS["analysis_skills"][t]["dependencies"]))
    skills = {}
    for name in [*targets, *deps]:
        path = find_skill(ws, name)
        skills[name] = {"installed": path is not None, "path": str(path) if path else None,
                        "kind": "analysis" if name in targets else "dependency"}
    licences = {}
    for name in deps:
        path = find_skill(ws, name)
        t = terms_of(path) if path else None
        if t:
            licences[name] = {"file": t[0], "links": t[1], "accepted": (ws / ".licenses" / t[0]).is_file()}
    keys = {}
    for t in targets:
        for k in TARGETS["analysis_skills"][t]["keys"]["required"]:
            keys[k] = {"required": True, "present": key_present(k)}
        for k in TARGETS["analysis_skills"][t]["keys"]["optional"]:
            keys.setdefault(k, {"required": False, "present": key_present(k)})
    uv = uv_path()
    npx = shutil.which("npx")
    project_env = (ws / ".env").is_file()
    return {
        "workspace": str(ws), "targets": targets,
        "system": {"os": platform.system(), "release": platform.release(), "version": platform.version(), **arch()},
        "uv": {"path": uv, "version": _version([uv, "--version"]) if uv else None},
        "node": {"node": _version(["node", "--version"]) if shutil.which("node") else None,
                 "npx": npx, "winget": bool(shutil.which("winget")), "brew": bool(shutil.which("brew"))},
        "agent_folders": {str(d): d.is_dir() for d in skill_dirs(ws)},
        "skills": skills, "licences": licences, "keys": keys,
        "project_env_file": project_env,
        "disk_free_gb": round(shutil.disk_usage(ws).free / 1e9, 1),
        "network": network(TARGETS["network"]) if net else {},
    }


def cmd_diagnose(args: argparse.Namespace) -> None:
    d = diagnose_data(args.ws, args.target, net=not args.no_network)
    if args.json:
        print(json.dumps(d, indent=2, ensure_ascii=False))
        return
    s = d["system"]
    print(f"System: {s['os']} {s['release']} ({s['native']}; Python seen by this script: {s['python_bits']} "
          f"{s['python_machine']})")
    print(f"uv: {d['uv']['version'] or 'NOT INSTALLED'}")
    print(f"Node.js / npx: {d['node']['node'] or 'NOT INSTALLED'}")
    print(f"Analysis skills prepared: {', '.join(d['targets'])}")
    for name, v in d["skills"].items():
        print(f"  {'ok     ' if v['installed'] else 'missing'} {v['kind']:<10} {name}")
    for name, v in d["licences"].items():
        print(f"  terms {'accepted' if v['accepted'] else 'NOT accepted'}: {name}")
    for k, v in d["keys"].items():
        print(f"  key {k}: {'set' if v['present'] else 'not set'}{'' if v['required'] else ' (optional)'}")
    if d["project_env_file"]:
        print("  note: there is a .env file inside the workspace; the analysis reads keys from ~/.env (home folder).")
    print(f"Free disk: {d['disk_free_gb']} GB")
    for u, r in d["network"].items():
        print(f"  network {u}: {r}")
    print("Next: setup.py plan (what is missing, as questions).")


# ---------------------------------------------------------------------------------------------- plan


def _item(key, title, question, why, identifiers, options, recommended, handled_by, agent_note=""):
    return {"key": key, "title": title, "question": question, "why": why, "identifiers": identifiers,
            "options": options, "recommended": recommended, "handled_by": handled_by, "agent_note": agent_note}


LATER = {"code": "later", "label": "Not now", "effect": "Nothing is installed; the analysis cannot run until this is in place."}


def plan_items(d: dict[str, Any], agents: list[str]) -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []
    sysname = d["system"]["os"]
    if not d["uv"]["path"]:
        items.append(_item(
            "setup.uv", "Install uv (runs the analysis scripts)",
            "uv is not installed. It runs every script of the analysis in its own environment with fixed package "
            "versions and brings its own Python, without changing the computer's Python. Install it for this user?",
            "Without uv no script of the analysis can run.",
            {"program": "uv (Astral)", "source": "https://docs.astral.sh/uv/", "administrator": "not needed"},
            [{"code": "install", "label": "Install uv", "effect": "The official installer puts uv in this user's folder "
              "(a few MB); the agent may need a new terminal session to see it."}, LATER], "install", "setup.py"))
    missing = [n for n, v in d["skills"].items() if not v["installed"]]
    if missing and not d["node"]["npx"]:
        how = ("winget (Windows package manager)" if sysname == "Windows" and d["node"]["winget"] else
               "Homebrew" if sysname == "Darwin" and d["node"]["brew"] else "the installer from nodejs.org")
        items.append(_item(
            "setup.node", "Install Node.js (installs the skills)",
            f"Node.js is not installed. Its `npx skills add` command downloads the {len(missing)} missing skill(s) "
            f"into this project. Install Node.js (LTS) with {how}?",
            "The skills are installed with npx, which comes with Node.js.",
            {"program": "Node.js LTS", "source": NODE_PAGE, "with": how},
            [{"code": "install", "label": "Install Node.js",
              "effect": ("It is installed by " + how + "; Windows may ask for permission.") if how != "the installer from nodejs.org"
              else "The download page opens; run the installer and tell the agent when it is done."}, LATER],
            "install", "setup.py"))
    for name in missing:
        v = d["skills"][name]
        if v["kind"] == "analysis":
            repo = TARGETS["analysis_skills"][name]["repo"]
            priv = TARGETS["analysis_skills"][name].get("private")
            items.append(_item(
                f"setup.skill:{name}", f"Install the analysis skill {name}",
                f"The analysis skill {name} is not installed. Install it from {repo}?",
                "It is the skill that runs the study.",
                {"skill": name, "repository": repo, "agents": ", ".join(agents)},
                [{"code": "install", "label": f"Install {name}",
                  "effect": f"npx skills add {repo} --skill {name} copies it into this project."}, LATER],
                "install", "npx",
                agent_note=("The repository is private: the researcher needs a GitHub account with access and to be "
                            "signed in (e.g. `gh auth login --web`, then `gh auth setup-git`) before installing.") if priv else ""))
        else:
            dep = TARGETS["dependency_skills"][name]
            items.append(_item(
                f"setup.skill:{name}", f"Install the skill {name}",
                f"The analysis uses the skill {name} ({dep['repo']}), which {dep['purpose']}. Install it in this project?",
                "Installing a skill brings third-party code and, for some, terms of use; each one is approved on its own.",
                {"skill": name, "repository": dep["repo"], "agents": ", ".join(agents)},
                [{"code": "install", "label": f"Install {name}",
                  "effect": f"npx skills add {dep['repo']} --skill {name} copies it into this project (latest version)."},
                 LATER], "install", "npx"))
    for name, v in d["licences"].items():
        if not v["accepted"]:
            items.append(_item(
                f"setup.licence:{name}", f"Terms of use of {name}",
                f"The skill {name} has terms of use that must be accepted before it is used. The agent shows them to you "
                "now. Do you accept them?",
                "The notice file records that these terms were seen and accepted; one acceptance never covers another skill.",
                {"skill": name, "terms": " ".join(v["links"]), "notice file": f".licenses/{v['file']}"},
                [{"code": "accept", "label": "I have read and accept the terms",
                  "effect": f"The agent records the acceptance in .licenses/{v['file']}."},
                 {"code": "decline", "label": "I do not accept",
                  "effect": f"Nothing is recorded; the steps that use {name} cannot run."}],
                None, "setup.py", agent_note=f"Show the terms first: setup.py terms {name}."))
    for k, v in d["keys"].items():
        if v["required"] and not v["present"]:
            items.append(_item(
                f"setup.key:{k}", "AlphaGenome access key" if k == "ALPHAGENOME_API_KEY" else f"Key {k}",
                f"The key {k} is not set. It is needed to fetch the predictions. Set it now? You paste it into a text "
                "file that opens on your screen; it never goes into the chat.",
                "The key identifies you to the service; it is kept in ~/.env and never shown.",
                {"variable": k, "stored in": "~/.env", "register": "https://deepmind.google.com/science/alphagenome/"},
                [{"code": "set", "label": "Set the key",
                  "effect": "~/.env opens in the text editor with an empty line for the key; paste it after '=', save, "
                            "close and tell the agent it is done."}, LATER],
                "set", "setup.py"))
    return items


def cmd_plan(args: argparse.Namespace) -> None:
    d = diagnose_data(args.ws, args.target, net=False)
    items = plan_items(d, args.agents)
    if args.json:
        print(json.dumps({"ready": not items, "items": items}, indent=2, ensure_ascii=False))
        sys.exit(1 if items else 0)
    if not items:
        print("Nothing missing. Next: setup.py verify.")
        return
    print(f"{len(items)} item(s) to go through with the researcher, one at a time, in their language. Never show the "
          "codes; install only after their approval of that item, with: setup.py install <key> --approved <key>.")
    for n, it in enumerate(items, 1):
        print(f"\n[{n}] {it['title']}\n    Question: {it['question']}\n    Why it matters: {it['why']}")
        print("    About: " + "; ".join(f"{k}: {v}" for k, v in it["identifiers"].items()))
        if it["agent_note"]:
            print(f"    For the agent: {it['agent_note']}")
        for o in it["options"]:
            rec = " (recommended)" if o["code"] == it["recommended"] else ""
            print(f"      - {o['label']}{rec}: {o['effect']}  [code: {o['code']}]")
        print(f"    Key: {it['key']}")
    sys.exit(1)


# --------------------------------------------------------------------------------------------- terms


def cmd_terms(args: argparse.Namespace) -> None:
    path = find_skill(args.ws, args.skill)
    if not path:
        sys.exit(f"The skill {args.skill} is not installed; install it first.")
    t = terms_of(path)
    if not t:
        print(f"The skill {args.skill} asks for no terms of use.")
        return
    print(f"Terms of use of {args.skill}: {' and '.join(t[1]) or '(see its SKILL.md)'}\n"
          f"Show these links to the researcher and ask whether they read and accept them. On their yes: "
          f"setup.py install setup.licence:{args.skill} --approved setup.licence:{args.skill}")


# ------------------------------------------------------------------------------------------- install


def _install_uv() -> int:
    cmd = UV_INSTALL["Windows" if platform.system() == "Windows" else "other"]
    print("Running the official uv installer: " + " ".join(cmd))
    code, out = _run(cmd, timeout=600)
    print(out[-2000:])
    found = uv_path()
    print(f"uv: {_version([found, '--version']) if found else 'not found yet'}"
          + ("" if shutil.which("uv") else f" (at {found}; a new terminal session puts it on PATH)" if found else ""))
    return 0 if found else (code or 1)


def _install_node() -> int:
    sysname = platform.system()
    if sysname == "Windows" and shutil.which("winget"):
        cmd = ["winget", "install", "-e", "--id", "OpenJS.NodeJS.LTS", "--accept-package-agreements",
               "--accept-source-agreements"]
    elif sysname == "Darwin" and shutil.which("brew"):
        cmd = ["brew", "install", "node"]
    else:
        webbrowser.open(NODE_PAGE)
        print(f"The Node.js download page was opened ({NODE_PAGE}). The researcher runs the installer (LTS) and tells "
              "the agent when it is done; then run setup.py diagnose again. On a server without a desktop, see "
              "references/platforms.md (user-level install with fnm).")
        return 0
    print("Running: " + " ".join(cmd))
    code, out = _run(cmd, timeout=1200)
    print(out[-2000:])
    print("Node.js: " + (_version(["node", "--version"]) or "not on PATH yet (open a new terminal session, or restart "
                                                              "the agent app, then run setup.py diagnose)"))
    return code


def _install_skill(ws: Path, name: str, agents: list[str]) -> int:
    if name in TARGETS["analysis_skills"]:
        repo = TARGETS["analysis_skills"][name]["repo"]
    elif name in TARGETS["dependency_skills"]:
        repo = TARGETS["dependency_skills"][name]["repo"]
    else:
        print(f"Unknown skill {name}.")
        return 2
    if find_skill(ws, name):
        print(f"{name} is already installed at {find_skill(ws, name)}; nothing was changed.")
        return 0
    npx = shutil.which("npx")
    if not npx:
        print("npx (Node.js) is not available; install Node.js first (setup.node).")
        return 2
    cmd = [npx, "--yes", "skills", "add", repo, "--skill", name, "--agent", *agents, "--copy", "-y", "--json"]
    print("Running: " + " ".join(["npx", *cmd[1:]]))
    code, out = _run(cmd, cwd=ws, timeout=900)
    m = re.search(r"\[\s*\{.*\}\s*\]", out, re.S)
    results = []
    if m:
        try:
            results = json.loads(m.group(0))
        except json.JSONDecodeError:
            results = []
    ok = [r for r in results if r.get("status") == "installed"]
    if not ok:
        print(out[-3000:])
        print(f"{name} was not installed (exit {code}).")
        return code or 1
    st = state(ws)
    for r in ok:
        st["installs"].append({"skill": r.get("name"), "source": r.get("source"), "hash": r.get("hash"),
                               "path": r.get("path"), "agents": r.get("agents"), "installed_at": _now()})
        print(f"Installed {r.get('name')} from {r.get('source')} at {r.get('path')}")
    save_state(ws, st)
    return 0


def _record_licence(ws: Path, name: str) -> int:
    path = find_skill(ws, name)
    t = terms_of(path) if path else None
    if not t:
        print(f"{name} is not installed or asks for no terms.")
        return 2
    target = ws / ".licenses" / t[0]
    target.parent.mkdir(exist_ok=True)
    target.write_text(f"Notice shown to the user on {_now()}: please check the terms of use of {name} at "
                      f"{' and '.join(t[1])}. The user confirmed they read and accept them.\n", encoding="utf-8")
    print(f"Recorded the acceptance in {target}")
    return 0


def _open_key_file(name: str) -> int:
    env = Path.home() / ".env"
    if key_present(name):
        print(f"{name} is already set; nothing was changed.")
        return 0
    lines = env.read_text(encoding="utf-8", errors="replace").splitlines() if env.is_file() else []
    if not any(line.strip().removeprefix("export ").strip().startswith(f"{name}=") for line in lines):
        with env.open("a", encoding="utf-8", newline="\n") as fh:  # append the empty line; existing lines are not shown
            if lines and lines[-1].strip():
                fh.write("\n")
            fh.write(f"{name}=\n")
    sysname = platform.system()
    if sysname == "Windows":
        subprocess.Popen(["notepad.exe", str(env)])
    elif sysname == "Darwin":
        subprocess.Popen(["open", "-e", str(env)])
    elif shutil.which("xdg-open"):
        subprocess.Popen(["xdg-open", str(env)])
    else:
        print(f"No text editor could be opened here. The researcher opens {env} (e.g. `nano ~/.env`), pastes the "
              f"key after '{name}=', saves and closes.")
        return 0
    print(f"{env} is open in the text editor. The researcher pastes the key right after '{name}=', saves and closes, "
          f"then says it is done; then run setup.py key {name}.")
    return 0


def cmd_install(args: argparse.Namespace) -> None:
    if args.approved != args.item:
        sys.exit(f"Not installed: --approved must name the same item ({args.item}) the researcher approved.")
    kind, _, subject = args.item.partition(":")
    if kind == "setup.uv":
        code = _install_uv()
    elif kind == "setup.node":
        code = _install_node()
    elif kind == "setup.skill":
        code = _install_skill(args.ws, subject, args.agents)
    elif kind == "setup.licence":
        code = _record_licence(args.ws, subject)
    elif kind == "setup.key":
        code = _open_key_file(subject)
    else:
        sys.exit(f"Unknown item {args.item}.")
    sys.exit(code)


def cmd_key(args: argparse.Namespace) -> None:
    present = key_present(args.name)
    print(f"{args.name}: {'set' if present else 'not set'}")
    sys.exit(0 if present else 1)


# ------------------------------------------------------------------------------------- verify, report


def cmd_verify(args: argparse.Namespace) -> None:
    uv = uv_path()
    if not uv:
        sys.exit("uv is not installed; run setup.py plan.")
    results = []
    for name in targets_for(args.ws, args.target):
        path = find_skill(args.ws, name)
        spec = TARGETS["analysis_skills"][name]
        if not path:
            results.append({"skill": name, "step": "installed", "ok": False, "detail": "not installed"})
            continue
        steps = [("check", spec["check"])]
        if (path / "assets" / "reference" / spec["reference"]).is_dir():
            steps.append(("reproduce", spec["reproduce"]))
        else:
            results.append({"skill": name, "step": "reproduce", "ok": None,
                            "detail": "no bundled reference study in this copy; skipped"})
        for step, argv in steps:
            t0 = time.monotonic()
            code, out = _run([uv, "run", "-q", str(path / argv[0]), *argv[1:]], cwd=args.ws, timeout=1800)
            tail = [l for l in out.strip().splitlines() if l.strip()][-3:]
            results.append({"skill": name, "step": step, "ok": code == 0, "seconds": round(time.monotonic() - t0, 1),
                            "detail": " | ".join(tail)[:600]})
    st = state(args.ws)
    st["verify"] = {"at": _now(), "results": results}
    save_state(args.ws, st)
    if args.json:
        print(json.dumps(results, indent=2, ensure_ascii=False))
    else:
        for r in results:
            mark = "ok  " if r["ok"] else ("skip" if r["ok"] is None else "FAIL")
            print(f"[{mark}] {r['skill']} {r['step']}" + (f" ({r['seconds']} s)" if "seconds" in r else "")
                  + f": {r['detail']}")
    sys.exit(0 if all(r["ok"] is not False for r in results) else 1)


def cmd_report(args: argparse.Namespace) -> None:
    d = diagnose_data(args.ws, args.target, net=False)
    st = state(args.ws)
    s = d["system"]
    lines = ["# Ambiente das skills AlphaGenome", "", f"Gerado em {_now()} por `alphagenome-skills-setup`.", "",
             "| Item | Situação |", "|---|---|",
             f"| Sistema | {s['os']} {s['release']} ({s['native']}) |",
             f"| uv | {d['uv']['version'] or 'não instalado'} |",
             f"| Node.js | {d['node']['node'] or 'não instalado'} |",
             f"| Espaço livre | {d['disk_free_gb']} GB |", "", "## Skills", "", "| Skill | Tipo | Instalada | Pasta |",
             "|---|---|---|---|"]
    lines += [f"| {n} | {'análise' if v['kind'] == 'analysis' else 'dependência'} | {'sim' if v['installed'] else 'não'} | "
              f"{v['path'] or '—'} |" for n, v in d["skills"].items()]
    if st["installs"]:
        lines += ["", "Instaladas por este setup (`npx skills add`; o hash identifica a versão baixada):", "",
                  "| Skill | Origem | Hash | Data |", "|---|---|---|---|"]
        lines += [f"| {i['skill']} | {i['source']} | `{(i['hash'] or '')[:16]}` | {i['installed_at']} |" for i in st["installs"]]
    lines += ["", "## Termos de uso e chaves", ""]
    lines += [f"- Termos de {n}: {'aceitos' if v['accepted'] else 'não aceitos'} ({' '.join(v['links'])})"
              for n, v in d["licences"].items()]
    lines += [f"- Chave {k}: {'definida' if v['present'] else 'não definida'}{'' if v['required'] else ' (opcional)'}"
              for k, v in d["keys"].items()]
    if st.get("verify"):
        lines += ["", f"## Verificação ({st['verify']['at']})", "", "| Skill | Etapa | Resultado | Detalhe |", "|---|---|---|---|"]
        lines += [f"| {r['skill']} | {r['step']} | {'ok' if r['ok'] else ('pulado' if r['ok'] is None else 'falhou')} | "
                  f"{r['detail'].replace('|', '/')} |" for r in st["verify"]["results"]]
    out = args.ws / REPORT_FILE
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


# -------------------------------------------------------------------------------------------------- cli


def main() -> None:
    for stream in (sys.stdout, sys.stderr):
        stream.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--workspace", default=".", help="The researcher's workspace (default: the current folder)")
    parser.add_argument("--target", action="append", help="Analysis skill to prepare (repeatable; default: those "
                                                          "installed, else expression-analysis)")
    parser.add_argument("--agents", nargs="+", default=["claude-code", "opencode"],
                        help="npx skills agent ids to install for (default: claude-code opencode, i.e. .claude/skills "
                             "and .agents/skills)")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("diagnose", help="What the computer has")
    p.add_argument("--json", action="store_true")
    p.add_argument("--no-network", action="store_true", help="Skip the network checks")
    p.set_defaults(func=cmd_diagnose)
    p = sub.add_parser("plan", help="What is missing, as questions")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_plan)
    p = sub.add_parser("terms", help="Show the terms of use a dependency skill asks for")
    p.add_argument("skill")
    p.set_defaults(func=cmd_terms)
    p = sub.add_parser("install", help="Install or record one item the researcher approved")
    p.add_argument("item", help="setup.uv, setup.node, setup.skill:<name>, setup.licence:<skill>, setup.key:<VAR>")
    p.add_argument("--approved", required=True, help="The same item, as approved by the researcher")
    p.set_defaults(func=cmd_install)
    p = sub.add_parser("key", help="Whether a key is set (never its value)")
    p.add_argument("name", nargs="?", default="ALPHAGENOME_API_KEY")
    p.set_defaults(func=cmd_key)
    p = sub.add_parser("verify", help="check (and offline reproduce) of each analysis skill")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_verify)
    p = sub.add_parser("report", help="Write environment-report.md")
    p.set_defaults(func=cmd_report)
    args = parser.parse_args()
    args.ws = Path(args.workspace).resolve()
    args.func(args)


if __name__ == "__main__":
    main()

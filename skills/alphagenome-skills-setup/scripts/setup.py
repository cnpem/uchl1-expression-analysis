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
  diagnose   What the computer has: system and architecture, uv, Node.js (npx), the agent app(s) and the
             folders they read skills from, the analysis and dependency skills (and copies kept where no
             chosen agent reads them), terms of use accepted, keys (presence only), disk space, network.
  plan       What is missing, as questions (title, question, why, options with label and effect,
             recommended, handled_by), in the shape of `expression.py check --json`.
  terms      Shows the terms of use a dependency skill asks to be accepted (links), for the agent to
             show before asking.
  install    Installs or records one approved item: setup.agent:<id> (records the agent app),
             setup.uv, setup.node, setup.skill:<name> (npx skills add), setup.place:<name> (copies a
             copy found elsewhere into the chosen agents' folders), setup.licence:<skill>, setup.key:<VAR> (opens the page where the key is created and the
             workspace's .env in the text editor, with an empty line for it and .env added to .gitignore;
             the researcher pastes the key there; it never reaches the chat).
  key        Says whether a key is set in the project's .env (the only place keys are read from), never
             its value.
  verify     Checks that each analysis skill (and this one) is in the folder of every chosen agent, runs
             its `check`, and its offline `reproduce` when the reference study is bundled; reports the
             services that cannot be reached and what they affect.
  report     Writes environment-report.md in the workspace: what is installed, versions, the result
             of verify; for the researcher to keep or attach when asking for help.

Run from the researcher's workspace (the folder the agent works in). Writes only: skills installed by
`npx skills add` or placed from a copy (into the chosen agents' skill folders of the workspace), .licenses/ files the researcher
accepted, the workspace's .env (only when asked, with an empty key line) and its .gitignore entry,
.alphagenome-setup/state.json and environment-report.md. Never prints, logs or copies the value of a key.
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
TERMS = json.loads((SKILL_DIR / "assets" / "terms.json").read_text(encoding="utf-8"))
KEY_PAGES = {"ALPHAGENOME_API_KEY": "https://alphagenome.google/api",
             "NCBI_API_KEY": "https://www.ncbi.nlm.nih.gov/account/settings/"}
AGENTS: dict[str, dict[str, str]] = TARGETS["agents"]
DEFAULT_AGENTS = ["claude-code", "opencode"]
SETUP_NAME = "alphagenome-skills-setup"
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


def _expand(folder: str) -> Path:
    """A folder of the agents table: `~` is the home folder; XDG_CONFIG_HOME and CLAUDE_CONFIG_DIR are honoured
    as npx skills does."""
    for prefix, var in (("~/.config/", "XDG_CONFIG_HOME"), ("~/.claude/", "CLAUDE_CONFIG_DIR")):
        if folder.startswith(prefix) and os.environ.get(var):
            return Path(os.environ[var]) / folder[len(prefix):]
    return Path(os.path.expanduser(folder))


def agent_dirs(ws: Path, agent: str) -> list[Path]:
    """The two folders an agent app reads skills from: the project's and the home folder's."""
    a = AGENTS[agent]
    return list(dict.fromkeys(p.resolve() for p in (ws / a["project"], _expand(a["global"]))))


def skill_dirs(ws: Path) -> list[Path]:
    """Every agent folder of the table (to find copies a chosen agent does not read)."""
    return list(dict.fromkeys(d for a in AGENTS for d in agent_dirs(ws, a)))


def locate(ws: Path, name: str, agents: list[str]) -> dict[str, Any]:
    """Where a skill is for each chosen agent, and the copies none of them reads.

    A skill counts as installed only when every chosen agent has it in one of its own folders: an agent lists
    the skills of its folders and no others, so a copy elsewhere (next to this skill, in `.skills/`, or in
    another agent's folder) is never loaded and the agent answers without it."""
    per = {a: next((str(d / name) for d in agent_dirs(ws, a) if (d / name / "SKILL.md").is_file()), None)
           for a in agents}
    seen = {v for v in per.values() if v}
    stray: list[str] = []
    for d in dict.fromkeys(c.resolve() for c in (SKILL_DIR.parent, ws / ".skills", ws / "skills", ws, *skill_dirs(ws))):
        f = str(d / name)
        if (d / name / "SKILL.md").is_file() and f not in seen and f not in stray:
            stray.append(f)
    return {"installed": all(per.values()), "path": next((v for v in per.values() if v), None),
            "per_agent": per, "stray": stray}


def find_skill(ws: Path, name: str, agents: list[str]) -> Path | None:
    """Any copy of a skill (in a chosen agent's folder first): enough to read its SKILL.md."""
    loc = locate(ws, name, agents)
    found = loc["path"] or (loc["stray"][0] if loc["stray"] else None)
    return Path(found) if found else None


def _agent_names(agents: list[str]) -> str:
    return ", ".join(AGENTS[a]["name"] for a in agents)


def _project_dirs(agents: list[str]) -> str:
    return " and ".join(dict.fromkeys(AGENTS[a]["project"] for a in agents))


def _key_in(path: Path, name: str) -> bool:
    try:
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip().removeprefix("export ").strip()
            key, sep, value = line.partition("=")
            if sep and key.strip() == name and value.strip().strip("\"'"):
                return True
    except OSError:
        pass
    return False


def key_where(ws: Path, name: str) -> str | None:
    """'project .env' when the project's .env sets the key, else None. The project's .env is the only place the
    analysis skills read keys from (not ~/.env, not environment variables). The value is never kept or shown."""
    return "project .env" if _key_in(ws / ".env", name) else None


def key_elsewhere(name: str) -> list[str]:
    """Places that hold the key but are not used (to tell the researcher), never their values."""
    out = []
    if os.environ.get(name):
        out.append("an environment variable")
    if _key_in(Path.home() / ".env", name):
        out.append("~/.env (home folder)")
    return out


def key_present(name: str, ws: Path | None = None) -> bool:
    return key_where(ws or Path.cwd(), name) is not None


def env_ignored(ws: Path) -> bool | None:
    """Whether git ignores the workspace's .env (None when the workspace is not a git repository)."""
    if not (ws / ".git").exists():
        return None
    code, _ = _run(["git", "check-ignore", "-q", ".env"], cwd=ws, timeout=30)
    return code == 0


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


def network(entries: list[dict[str, str]]) -> dict[str, dict[str, Any]]:
    """Whether each service answers (any HTTP status counts as reachable), with what it affects."""
    out = {}
    for e in entries:
        try:
            req = urllib.request.Request(e["url"], headers={"User-Agent": "alphagenome-skills-setup"})
            with urllib.request.urlopen(req, timeout=10) as r:
                status, ok = f"ok ({r.status})", True
        except urllib.error.HTTPError as exc:  # the server answered: reachable
            status, ok = f"ok ({exc.code})", True
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            status, ok = f"unreachable ({getattr(exc, 'reason', exc)})", False
        out[e["url"]] = {"service": e["service"], "affects": e["affects"], "status": status, "reachable": ok}
    return out


def unreachable_notes(net: dict[str, dict[str, Any]]) -> list[str]:
    return [f"{v['service']} cannot be reached now ({v['status']}): this affects {v['affects']}. It is usually the "
            "service, not this computer; try again later." for v in net.values() if not v["reachable"]]


def state(ws: Path) -> dict[str, Any]:
    p = ws / STATE_DIR / "state.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else {"installs": [], "verify": None}


def chosen_agents(ws: Path, given: list[str] | None) -> tuple[list[str], str]:
    """The agent apps to prepare and where that came from: --agents, the record of an earlier run, the agents
    of earlier installs, or the default (nobody asked yet)."""
    if given:
        bad = [a for a in given if a not in AGENTS]
        if bad:
            sys.exit(f"Unknown agent id(s) {bad}; known: {sorted(AGENTS)}")
        return list(dict.fromkeys(given)), "given"
    st = state(ws)
    if st.get("agents"):
        return list(st["agents"]), "recorded"
    names = {n for i in st.get("installs", []) for n in (i.get("agents") or [])}
    inferred = [a for a, v in AGENTS.items() if any(n and n in v["name"] for n in names)]
    if inferred:
        return inferred, "earlier installs"
    return list(DEFAULT_AGENTS), "default"


def record_agents(ws: Path, agents: list[str]) -> None:
    st = state(ws)
    st["agents"] = list(dict.fromkeys([*(st.get("agents") or []), *agents]))
    save_state(ws, st)


def save_state(ws: Path, st: dict[str, Any]) -> None:
    (ws / STATE_DIR).mkdir(exist_ok=True)
    (ws / STATE_DIR / "state.json").write_text(json.dumps(st, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def targets_for(ws: Path, wanted: list[str] | None, agents: list[str]) -> list[str]:
    """Analysis skills to prepare: the ones named, else the ones installed, else expression-analysis."""
    names = list(TARGETS["analysis_skills"])
    if wanted:
        bad = [w for w in wanted if w not in names]
        if bad:
            sys.exit(f"Unknown analysis skill(s) {bad}; known: {names}")
        return wanted
    def here(n: str) -> bool:  # in a chosen agent's folder, in this project, or next to this skill
        loc = locate(ws, n, agents)
        return bool(loc["path"]) or any(Path(s).is_relative_to(ws) or Path(s).parent == SKILL_DIR.parent.resolve()
                                        for s in loc["stray"])
    return [n for n in names if here(n)] or ["expression-analysis"]


# ------------------------------------------------------------------------------------------ diagnose


def diagnose_data(ws: Path, wanted: list[str] | None, agents: list[str], agents_source: str,
                  net: bool = True) -> dict[str, Any]:
    targets = targets_for(ws, wanted, agents)
    deps = list(dict.fromkeys(d for t in targets for d in TARGETS["analysis_skills"][t]["dependencies"]))
    skills = {}
    for name in [*targets, *deps]:
        skills[name] = {**locate(ws, name, agents), "kind": "analysis" if name in targets else "dependency"}
    licences = {}
    for name in deps:
        path = find_skill(ws, name, agents)
        t = terms_of(path) if path else None
        if t:
            licences[name] = {"file": t[0], "links": t[1], "accepted": (ws / ".licenses" / t[0]).is_file()}
    keys = {}
    for t in targets:
        for k in TARGETS["analysis_skills"][t]["keys"]["required"]:
            keys[k] = {"required": True, "present": key_present(k, ws), "where": key_where(ws, k),
                       "unused_elsewhere": key_elsewhere(k)}
        for k in TARGETS["analysis_skills"][t]["keys"]["optional"]:
            keys.setdefault(k, {"required": False, "present": key_present(k, ws), "where": key_where(ws, k),
                                "unused_elsewhere": key_elsewhere(k)})
    uv = uv_path()
    npx = shutil.which("npx")
    project_env = (ws / ".env").is_file()
    return {
        "workspace": str(ws), "targets": targets,
        "agents": {"ids": agents, "source": agents_source,
                   "folders": {a: [str(d) for d in agent_dirs(ws, a)] for a in agents}},
        "setup_skill": {"path": str(SKILL_DIR), **locate(ws, SETUP_NAME, agents)},
        "system": {"os": platform.system(), "release": platform.release(), "version": platform.version(), **arch()},
        "uv": {"path": uv, "version": _version([uv, "--version"]) if uv else None},
        "node": {"node": _version(["node", "--version"]) if shutil.which("node") else None,
                 "npx": npx, "winget": bool(shutil.which("winget")), "brew": bool(shutil.which("brew"))},
        "skills": skills, "licences": licences, "keys": keys,
        "project_env_file": project_env,
        "project_env_ignored_by_git": env_ignored(ws) if project_env else None,
        "disk_free_gb": round(shutil.disk_usage(ws).free / 1e9, 1),
        "network": network(TARGETS["network"]) if net else {},
    }


def cmd_diagnose(args: argparse.Namespace) -> None:
    d = diagnose_data(args.ws, args.target, args.agents, args.agents_source, net=not args.no_network)
    if args.json:
        print(json.dumps(d, indent=2, ensure_ascii=False))
        return
    s = d["system"]
    print(f"System: {s['os']} {s['release']} ({s['native']}; Python seen by this script: {s['python_bits']} "
          f"{s['python_machine']})")
    print(f"uv: {d['uv']['version'] or 'NOT INSTALLED'}")
    print(f"Node.js / npx: {d['node']['node'] or 'NOT INSTALLED'}")
    ag = d["agents"]
    src = {"given": "chosen", "recorded": "recorded by an earlier run", "earlier installs": "from earlier installs",
           "default": "NOT ASKED YET: ask which app the researcher uses and pass --agents; until then"}[ag["source"]]
    print(f"Agent app(s) ({src}): {_agent_names(ag['ids'])}")
    for a, dirs in ag["folders"].items():
        print(f"  {AGENTS[a]['name']} reads skills from: {' or '.join(dirs)}")
    print(f"Analysis skills prepared: {', '.join(d['targets'])}")
    for name, v in d["skills"].items():
        print(f"  {'ok     ' if v['installed'] else 'missing'} {v['kind']:<10} {name}"
              + "".join(f"\n          not for {AGENTS[a]['name']}" for a, p in v["per_agent"].items() if not p and v["path"]))
        for s in v["stray"]:
            print(f"          a copy is at {s}, a folder the chosen agent app(s) do not read")
    if not d["setup_skill"]["installed"]:
        print(f"  WARNING: this setup skill is at {SKILL_DIR}, which {_agent_names(ag['ids'])} do not read: a new "
              "conversation will not list it. `plan` offers to place it.")
    for name, v in d["licences"].items():
        print(f"  terms {'accepted' if v['accepted'] else 'NOT accepted'}: {name}")
    for k, v in d["keys"].items():
        print(f"  key {k}: {('set (' + v['where'] + ')') if v['present'] else 'not set in the project .env'}"
              f"{'' if v['required'] else ' (optional)'}"
              + (f"; also found in {', '.join(v['unused_elsewhere'])}, which the analysis does not use"
                 if v["unused_elsewhere"] else ""))
    if d["project_env_file"] and d["project_env_ignored_by_git"] is False:
        print("  WARNING: the workspace's .env is not ignored by git; a commit would publish the key. "
              "`install setup.key:…` adds it to .gitignore.")
    print(f"Free disk: {d['disk_free_gb']} GB")
    for u, r in d["network"].items():
        print(f"  network {r['service']} ({u}): {r['status']}" + ("" if r["reachable"] else f"; affects {r['affects']}"))
    print("Next: setup.py plan (what is missing, as questions).")


# ---------------------------------------------------------------------------------------------- plan


def _item(key, title, question, why, identifiers, options, recommended, handled_by, agent_note=""):
    return {"key": key, "title": title, "question": question, "why": why, "identifiers": identifiers,
            "options": options, "recommended": recommended, "handled_by": handled_by, "agent_note": agent_note}


LATER = {"code": "later", "label": "Not now", "effect": "Nothing is installed; the analysis cannot run until this is in place."}


def _place_item(name: str, v: dict[str, Any], agents: list[str], kind: str) -> dict[str, Any]:
    src = v["path"] or v["stray"][0]
    missing = [a for a, p in v["per_agent"].items() if not p]
    return _item(
        f"setup.place:{name}", f"Put the {kind} {name} where your agent app reads it",
        f"A copy of {name} is at {src}, but {_agent_names(missing)} only reads skills from "
        f"{_project_dirs(missing)} (or its folder in your home). In a new conversation it would not know this skill "
        f"and would answer without it. Copy it into {_project_dirs(missing)}?",
        "An agent lists the skills of its own folders and no others; this is how it knows the skill exists.",
        {"skill": name, "copy found": src, "copied to": _project_dirs(missing), "agents": _agent_names(missing)},
        [{"code": "place", "label": "Copy it there", "effect": f"The folder is copied into {_project_dirs(missing)}; "
          f"the copy at {src} is left as it is."}, LATER], "place", "setup.py",
        agent_note=(f"To install it from the repository instead (if the researcher prefers), the item is "
                    f"setup.skill:{name}." if name in TARGETS["analysis_skills"] or name in TARGETS["dependency_skills"]
                    else ""))


def plan_items(d: dict[str, Any], agents: list[str]) -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []
    sysname = d["system"]["os"]
    if d["agents"]["source"] == "default":
        items.append(_item(
            "setup.agent", "Which app do you use with the skills?",
            "Each app reads skills from its own folder, so the skills go where your app looks. Which app do you use "
            "(you can name more than one)?",
            "A skill in another app's folder is not listed in yours, and the agent then works without it.",
            {"apps": ", ".join(v["name"] for v in AGENTS.values())},
            [{"code": a, "label": v["name"], "effect": f"Skills go to {v['project']} in this project."}
             for a, v in AGENTS.items()], None, "setup.py",
            agent_note="On the answer: setup.py install setup.agent:<id> --approved setup.agent:<id> (once per app), "
                       "then run plan again; nothing else is installed before this."))
        return items
    if not d["setup_skill"]["installed"]:
        items.append(_place_item(SETUP_NAME, d["setup_skill"], agents, "setup skill"))
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
    # an analysis skill found elsewhere is copied (it may be a version not yet published); dependency skills always
    # come from npx skills add, which records their source and hash
    placeable = [n for n in missing if d["skills"][n]["kind"] == "analysis"
                 and (d["skills"][n]["path"] or d["skills"][n]["stray"])]
    for name in placeable:
        items.append(_place_item(name, d["skills"][name], agents, f"{d['skills'][name]['kind']} skill"))
    missing = [n for n in missing if n not in placeable]
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
            info = TERMS.get(name, {})
            items.append(_item(
                f"setup.licence:{name}", f"Terms of use of {name}",
                f"The skill {name} uses {info.get('service', 'an external service')}, whose terms of use must be "
                "accepted before it is used. The agent explains the main points and opens the pages for you to read. "
                "Do you accept them?",
                "The notice file records that these terms were seen and accepted; one acceptance never covers another skill.",
                {"skill": name, "terms": " ".join(v["links"]), "notice file": f".licenses/{v['file']}"},
                [{"code": "accept", "label": "I have read and accept the terms",
                  "effect": f"The agent records the acceptance in .licenses/{v['file']}."},
                 {"code": "decline", "label": "I do not accept",
                  "effect": f"Nothing is recorded; the steps that use {name} cannot run."}],
                None, "setup.py",
                agent_note=f"Before asking: `setup.py terms {name} --open` (explain its points in the researcher's "
                           "language; the pages open in their browser) and give them time to read."))
    for k, v in d["keys"].items():
        if v["required"] and not v["present"]:
            items.append(_item(
                f"setup.key:{k}", "AlphaGenome access key" if k == "ALPHAGENOME_API_KEY" else f"Key {k}",
                f"The key {k} is not set. It is needed to fetch the predictions. Create it now? The agent guides you "
                "step by step: a page opens for you to sign in and create the key, and a text file opens for you to "
                "paste it; it never goes into the chat.",
                "The key identifies you to the service; it is personal, kept only in this project's .env (which git "
                "ignores) and never shown." + (f" A key in {', '.join(v['unused_elsewhere'])} is not used." if v["unused_elsewhere"] else ""),
                {"variable": k, "stored in": ".env of this project", "key page": KEY_PAGES.get(k, "—")},
                [{"code": "set", "label": "Create and set the key",
                  "effect": "The key page opens in the browser and .env in the text editor; .env is added to .gitignore."},
                 LATER],
                "set", "setup.py",
                agent_note="Guide the researcher with references/keys.md, one step at a time (account, sign in, terms, "
                           "create, copy, paste, save), before and after running the install."))
    return items


def cmd_plan(args: argparse.Namespace) -> None:
    d = diagnose_data(args.ws, args.target, args.agents, args.agents_source, net=not args.no_network)
    items = plan_items(d, args.agents)
    notes = unreachable_notes(d["network"])
    if args.json:
        print(json.dumps({"ready": not items, "agents": d["agents"], "notices": notes, "items": items},
                         indent=2, ensure_ascii=False))
        sys.exit(1 if items else 0)
    for n in notes:
        print(f"Notice (tell the researcher; nothing to install): {n}")
    if not items:
        print("Nothing missing. Next: setup.py verify, then setup.py report. The setup is finished only after both: "
              "verify is what shows that each skill is where the agent app reads it and that its check passes.")
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
    print("\nAfter the last item: setup.py plan again, then verify and report (the setup is not finished before).")
    sys.exit(1)


# --------------------------------------------------------------------------------------------- terms


def cmd_terms(args: argparse.Namespace) -> None:
    path = find_skill(args.ws, args.skill, args.agents)
    if not path:
        sys.exit(f"The skill {args.skill} is not installed; install it first.")
    t = terms_of(path)
    if not t:
        print(f"The skill {args.skill} asks for no terms of use.")
        return
    info = TERMS.get(args.skill, {})
    links = dict(info.get("links") or {})
    for u in t[1]:  # the links the skill itself names come first in its notice
        if u not in links.values():
            links[f"Named by {args.skill}"] = u
    print(f"Terms of use of {args.skill}" + (f": {info['service']}" if info.get("service") else ""))
    for point in info.get("points", []):
        print(f"  - {point}")
    for label, url in links.items():
        print(f"  {label}: {url}")
    if args.open:
        for url in links.values():
            webbrowser.open(url)
        print("The pages were opened in the browser.")
    print(f"Explain the points above in the researcher's language (they are a guide; the pages are what counts), "
          f"let them read, and ask whether they accept. On their yes: setup.py install setup.licence:{args.skill} "
          f"--approved setup.licence:{args.skill}")


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
    loc = locate(ws, name, agents)
    missing = [a for a, p in loc["per_agent"].items() if not p]
    if not missing:
        print(f"{name} is already in the folder of {_agent_names(agents)} ({loc['path']}); nothing was changed.")
        return 0
    npx = shutil.which("npx")
    if not npx:
        print("npx (Node.js) is not available; install Node.js first (setup.node).")
        return 2
    cmd = [npx, "--yes", "skills", "add", repo, "--skill", name, "--agent", *missing, "--copy", "-y", "--json"]
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
    record_agents(ws, agents)
    still = [a for a, p in locate(ws, name, agents)["per_agent"].items() if not p]
    if still:
        print(f"WARNING: npx reported the install, but {name} is not in the folder of {_agent_names(still)} "
              f"({_project_dirs(still)}); run diagnose.")
        return 1
    return 0


def _place_skill(ws: Path, name: str, agents: list[str]) -> int:
    """Copies a copy of the skill found elsewhere into the project folder of each chosen agent that lacks it."""
    loc = locate(ws, name, agents)
    src = loc["path"] or (loc["stray"][0] if loc["stray"] else None)
    if src is None:
        print(f"No copy of {name} was found; install it with setup.skill:{name}.")
        return 2
    missing = [a for a, p in loc["per_agent"].items() if not p]
    if not missing:
        print(f"{name} is already in the folder of {_agent_names(agents)}; nothing was changed.")
        return 0
    st = state(ws)
    for dest_dir in dict.fromkeys(ws / AGENTS[a]["project"] for a in missing):
        dest = dest_dir / name
        if (dest / "SKILL.md").is_file():
            continue
        if dest.exists():
            print(f"{dest} exists without a SKILL.md; it was left as it is. Remove or rename it, then run this again.")
            return 1
        shutil.copytree(src, dest, ignore=shutil.ignore_patterns("__pycache__", "__marimo__", "*.pyc", "*:Zone.Identifier"))
        print(f"Copied {src} to {dest}")
        st["installs"].append({"skill": name, "source": f"copy of {src}", "hash": None, "path": str(dest),
                               "agents": [AGENTS[a]["name"] for a in missing if ws / AGENTS[a]["project"] == dest_dir],
                               "installed_at": _now()})
    save_state(ws, st)
    record_agents(ws, agents)
    print(f"The copy at {src} was left as it is; the researcher may delete it. A new conversation of the agent app "
          "lists the skill.")
    return 0


def _record_licence(ws: Path, name: str, agents: list[str]) -> int:
    path = find_skill(ws, name, agents)
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


def _ensure_ignored(ws: Path) -> str:
    """Adds `.env` to the workspace's .gitignore (created if missing) unless git already ignores it."""
    if env_ignored(ws):
        return "already ignored by git"
    gi = ws / ".gitignore"
    lines = gi.read_text(encoding="utf-8", errors="replace").splitlines() if gi.is_file() else []
    if ".env" in (l.strip() for l in lines):
        return "listed in .gitignore"
    with gi.open("a", encoding="utf-8", newline="\n") as fh:
        if lines and lines[-1].strip():
            fh.write("\n")
        fh.write("# personal access keys (never commit)\n.env\n")
    return "added to .gitignore"


def _open_key_file(ws: Path, name: str) -> int:
    where = key_where(ws, name)
    if where:
        print(f"{name} is already set ({where}); nothing was changed.")
        return 0
    env = ws / ".env"
    ignored = _ensure_ignored(ws)
    lines = env.read_text(encoding="utf-8", errors="replace").splitlines() if env.is_file() else []
    if not any(line.strip().removeprefix("export ").strip().startswith(f"{name}=") for line in lines):
        with env.open("a", encoding="utf-8", newline="\n") as fh:  # append the empty line; existing lines are not shown
            if lines and lines[-1].strip():
                fh.write("\n")
            fh.write(f"{name}=\n")
    page = KEY_PAGES.get(name)
    if page:
        webbrowser.open(page)
    sysname = platform.system()
    if sysname == "Windows":
        subprocess.Popen(["notepad.exe", str(env)])
    elif sysname == "Darwin":
        subprocess.Popen(["open", "-e", str(env)])
    elif shutil.which("xdg-open"):
        subprocess.Popen(["xdg-open", str(env)])
    else:
        print(f"No text editor could be opened here. The researcher opens {env} (e.g. `nano .env`), pastes the key "
              f"after '{name}=', saves and closes.")
        return 0
    print((f"The key page ({page}) is open in the browser. " if page else "")
          + f"{env} is open in the text editor (.env {ignored}). The researcher creates the key, pastes it right "
            f"after '{name}=', saves and closes, then says it is done; then run setup.py key {name}.")
    return 0


def cmd_install(args: argparse.Namespace) -> None:
    if args.approved != args.item:
        sys.exit(f"Not installed: --approved must name the same item ({args.item}) the researcher approved.")
    kind, _, subject = args.item.partition(":")
    if kind == "setup.agent":
        if subject not in AGENTS:
            sys.exit(f"Unknown agent id {subject!r}; known: {sorted(AGENTS)}")
        record_agents(args.ws, [subject])
        print(f"Recorded: the skills go to {AGENTS[subject]['project']} for {AGENTS[subject]['name']}. Next: setup.py plan.")
        code = 0
    elif kind == "setup.uv":
        code = _install_uv()
    elif kind == "setup.node":
        code = _install_node()
    elif kind == "setup.skill":
        code = _install_skill(args.ws, subject, args.agents)
    elif kind == "setup.place":
        code = _place_skill(args.ws, subject, args.agents)
    elif kind == "setup.licence":
        code = _record_licence(args.ws, subject, args.agents)
    elif kind == "setup.key":
        code = _open_key_file(args.ws, subject)
    else:
        sys.exit(f"Unknown item {args.item}.")
    sys.exit(code)


def cmd_key(args: argparse.Namespace) -> None:
    where = key_where(args.ws, args.name)
    print(f"{args.name}: {'set (' + where + ')' if where else 'not set'}")
    if (args.ws / ".env").is_file() and env_ignored(args.ws) is False:
        print("WARNING: the workspace's .env is not ignored by git; add '.env' to .gitignore before any commit.")
    sys.exit(0 if where else 1)


# ------------------------------------------------------------------------------------- verify, report


def cmd_verify(args: argparse.Namespace) -> None:
    uv = uv_path()
    if not uv:
        sys.exit("uv is not installed; run setup.py plan.")
    results = []
    agents = args.agents
    if args.agents_source == "default":
        results.append({"skill": "-", "step": "agent app", "ok": False,
                        "detail": "not asked: run plan (its first item asks which app the researcher uses)"})
    me = locate(args.ws, SETUP_NAME, agents)
    if not me["installed"]:
        results.append({"skill": SETUP_NAME, "step": "placed", "ok": None,
                        "detail": f"at {SKILL_DIR}, not read by {_agent_names([a for a, p in me['per_agent'].items() if not p])}; "
                                  f"plan offers setup.place:{SETUP_NAME}"})
    for name in targets_for(args.ws, args.target, agents):
        loc = locate(args.ws, name, agents)
        spec = TARGETS["analysis_skills"][name]
        if not loc["installed"]:
            lacking = [a for a, p in loc["per_agent"].items() if not p]
            where = f"; a copy is at {loc['stray'][0]}" if loc["stray"] else ""
            results.append({"skill": name, "step": "placed", "ok": False,
                            "detail": f"not in {_project_dirs(lacking)}, the folder {_agent_names(lacking)} reads{where}; "
                                      "the agent app would not list it (plan offers the fix)"})
            continue
        path = Path(loc["path"])
        results.append({"skill": name, "step": "placed", "ok": True,
                        "detail": "; ".join(f"{AGENTS[a]['name']}: {p}" for a, p in loc["per_agent"].items())})
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
    if not args.no_network:
        for v in network(TARGETS["network"]).values():
            if not v["reachable"]:
                results.append({"skill": "-", "step": "network", "ok": None,
                                "detail": f"{v['service']} cannot be reached ({v['status']}); affects {v['affects']}"})
    st = state(args.ws)
    st["verify"] = {"at": _now(), "results": results}
    save_state(args.ws, st)
    if args.json:
        print(json.dumps(results, indent=2, ensure_ascii=False))
    else:
        for r in results:
            mark = "ok  " if r["ok"] else ("note" if r["ok"] is None else "FAIL")
            print(f"[{mark}] {r['skill']} {r['step']}" + (f" ({r['seconds']} s)" if "seconds" in r else "")
                  + f": {r['detail']}")
    sys.exit(0 if all(r["ok"] is not False for r in results) else 1)


def cmd_report(args: argparse.Namespace) -> None:
    d = diagnose_data(args.ws, args.target, args.agents, args.agents_source, net=False)
    st = state(args.ws)
    s = d["system"]
    lines = ["# Ambiente das skills AlphaGenome", "", f"Gerado em {_now()} por `alphagenome-skills-setup`.", "",
             "| Item | Situação |", "|---|---|",
             f"| Sistema | {s['os']} {s['release']} ({s['native']}) |",
             f"| uv | {d['uv']['version'] or 'não instalado'} |",
             f"| Node.js | {d['node']['node'] or 'não instalado'} |",
             f"| Espaço livre | {d['disk_free_gb']} GB |",
             f"| App(s) do agente | {_agent_names(d['agents']['ids'])}"
             + (" (não perguntado: padrão)" if d['agents']['source'] == 'default' else "") + " |",
             f"| Pastas de skills | {_project_dirs(d['agents']['ids'])} |",
             "", "## Skills", "", "| Skill | Tipo | Na pasta do agente | Pasta | Cópias fora das pastas do agente |",
             "|---|---|---|---|---|"]
    rows = [(SETUP_NAME, "setup", d["setup_skill"]), *((n, 'análise' if v['kind'] == 'analysis' else 'dependência', v)
                                                      for n, v in d["skills"].items())]
    lines += [f"| {n} | {k} | {'sim' if v['installed'] else 'NÃO'} | {v['path'] or '—'} | {', '.join(v['stray']) or '—'} |"
              for n, k, v in rows]
    if st["installs"]:
        lines += ["", "Instaladas por este setup (`npx skills add`; o hash identifica a versão baixada):", "",
                  "| Skill | Origem | Hash | Data |", "|---|---|---|---|"]
        lines += [f"| {i['skill']} | {i['source']} | `{(i['hash'] or '')[:16]}` | {i['installed_at']} |" for i in st["installs"]]
    lines += ["", "## Termos de uso e chaves", ""]
    lines += [f"- Termos de {n}: {'aceitos' if v['accepted'] else 'não aceitos'} ({' '.join(v['links'])})"
              for n, v in d["licences"].items()]
    lines += [f"- Chave {k}: {('definida (' + v['where'] + ')') if v['present'] else 'não definida'}"
              f"{'' if v['required'] else ' (opcional)'}" for k, v in d["keys"].items()]
    if d["project_env_file"]:
        lines.append(f"- `.env` do projeto ignorado pelo git: "
                     f"{ {True: 'sim', False: 'NÃO (risco de publicar a chave)', None: 'não é um repositório git'}[d['project_env_ignored_by_git']] }")
    if st.get("verify"):
        lines += ["", f"## Verificação ({st['verify']['at']})", "", "| Skill | Etapa | Resultado | Detalhe |", "|---|---|---|---|"]
        lines += [f"| {r['skill']} | {r['step']} | {'ok' if r['ok'] else ('aviso' if r['ok'] is None else 'falhou')} | "
                  f"{r['detail'].replace('|', '/')} |" for r in st["verify"]["results"]]
    else:
        lines += ["", "## Verificação", "", "Não rodou: o ambiente ainda não foi verificado (`setup.py verify`)."]
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
    parser.add_argument("--agents", action="append", default=None,
                        help="agent app id(s), comma-separated or repeated (assets/targets.json `agents`), e.g. "
                             "--agents antigravity or --agents claude-code,opencode: claude-code reads .claude/skills, the others "
                             ".agents/skills. Default: those recorded by an earlier run, else claude-code opencode "
                             "(and plan first asks which app the researcher uses)")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("diagnose", help="What the computer has")
    p.add_argument("--json", action="store_true")
    p.add_argument("--no-network", action="store_true", help="Skip the network checks")
    p.set_defaults(func=cmd_diagnose)
    p = sub.add_parser("plan", help="What is missing, as questions")
    p.add_argument("--json", action="store_true")
    p.add_argument("--no-network", action="store_true", help="Skip the network checks")
    p.set_defaults(func=cmd_plan)
    p = sub.add_parser("terms", help="Explain the terms of use a dependency skill asks for, with their links")
    p.add_argument("skill")
    p.add_argument("--open", action="store_true", help="Open the pages in the browser")
    p.set_defaults(func=cmd_terms)
    p = sub.add_parser("install", help="Install or record one item the researcher approved")
    p.add_argument("item", help="setup.uv, setup.node, setup.skill:<name>, setup.licence:<skill>, setup.key:<VAR>")
    p.add_argument("--approved", required=True, help="The same item, as approved by the researcher")
    p.set_defaults(func=cmd_install)
    p = sub.add_parser("key", help="Whether a key is set (never its value)")
    p.add_argument("name", nargs="?", default="ALPHAGENOME_API_KEY")
    p.set_defaults(func=cmd_key)
    p = sub.add_parser("verify", help="placement, check (and offline reproduce) of each analysis skill")
    p.add_argument("--json", action="store_true")
    p.add_argument("--no-network", action="store_true", help="Skip the network checks")
    p.set_defaults(func=cmd_verify)
    p = sub.add_parser("report", help="Write environment-report.md")
    p.set_defaults(func=cmd_report)
    args = parser.parse_args()
    args.ws = Path(args.workspace).resolve()
    given = [a.strip() for v in (args.agents or []) for a in v.split(",") if a.strip()]
    args.agents, args.agents_source = chosen_agents(args.ws, given or None)
    args.func(args)


if __name__ == "__main__":
    main()

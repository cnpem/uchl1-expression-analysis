---
name: alphagenome-skills-setup
description: "Prepares a researcher's computer (Windows, macOS or Linux) to use the AlphaGenome analysis skills expression-analysis and splicing-analysis: checks what is there, then, with the researcher's approval for each item, installs uv and Node.js, installs the analysis and dependency skills with npx skills add, records the terms of use the researcher accepts, opens ~/.env for them to paste the AlphaGenome key without it reaching the chat, verifies with each skill's check (and its offline reproduce when the reference is bundled) and writes an environment report. Use it whenever someone wants to install, set up or prepare the environment for expression-analysis or splicing-analysis, when one of their checks reports missing skills, terms or keys, or when a researcher without terminal experience asks how to start using these skills."
metadata:
  version: 0.1.0
---

# alphagenome-skills-setup

## Overview

The analysis skills need a few things on the computer: `uv` (runs their scripts with fixed package
versions and its own Python), their dependency skills (DeepMind science skills and the K-Dense
`alphagenome` skill), the terms of use of some of them accepted, and an AlphaGenome key. This skill
finds what is missing and installs it **one item at a time, only after the researcher approves that
item**, so a researcher who never used a terminal can get to a working setup by answering questions.

Why each rule matters:
- **One approval per item.** Installing brings third-party code, and some skills come with terms of
  use; the researcher decides each one, and the record (`.licenses/`, the report) shows what they agreed to.
- **The key never passes through the chat.** It is pasted by the researcher into `~/.env` in a text
  editor; scripts only check that it is there. Do not read, print, copy or `cat` `~/.env`.
- **No administrator rights, no Docker.** Everything installs for the user; Windows may still ask for
  permission when Node.js is installed with winget.
- **The agent's shell cannot type into prompts**, so every command here runs without questions.

`SETUP` stands for this skill's folder. Run every command from the researcher's workspace (the folder
the agent works in, where the study will live). Every command works unchanged in PowerShell, bash and zsh.

## Workflow

Talk in the researcher's language, in short sentences. Do not show codes or keys of the questions.

1. **Say what will happen** in two sentences: you will check the computer, then ask before installing
   each missing piece; it takes a few minutes, mostly downloads.
2. **uv first.** Run `uv --version`. If it fails, uv is missing and `setup.py` cannot run yet: explain
   that uv runs the analysis scripts with their own Python, ask to install it, and on their yes run the
   official installer for the system (`references/platforms.md`):
   - Windows (PowerShell): `powershell -NoProfile -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"`
   - macOS / Linux: `curl -LsSf https://astral.sh/uv/install.sh | sh`
   If `uv` is still not found, use its full path (`~/.local/bin/uv`, `%USERPROFILE%\.local\bin\uv.exe`).
3. **Which agent.** Ask which app they use (Claude, OpenCode, Antigravity, Codex, Cursor, Gemini CLI)
   and pass its id with `--agents` (the `agents` table of `assets/targets.json`); without an answer the
   default installs for Claude and OpenCode (`.claude/skills` and `.agents/skills`).
4. **Diagnose**: `uv run SETUP/scripts/setup.py diagnose`. Show a short table: what is there, what is
   missing, anything unusual (no network to a domain, little disk, a `.env` inside the project, Windows
   on ARM). Nothing is installed here.
5. **Plan**: `uv run SETUP/scripts/setup.py plan` (or `--json`). For each item, one at a time: its
   title, the question, why it matters, what it is about, each option's label and what happens, and
   the recommended one. Follow its "For the agent" line first (e.g. show the terms with `setup.py terms`).
6. **Install what they approve**: `uv run SETUP/scripts/setup.py install <key> --approved <key>`, with
   the item's key; one call per approved item. A "not now" is fine: say what will not work without it
   and go on. A failure does not stop the rest: explain it in plain words (`references/troubleshooting.md`)
   and offer the alternative.
   - After Node.js, a new terminal session (or restarting the agent app) may be needed before `npx` works.
   - After the dependency skills are installed, run `plan` again: their terms of use only appear then.
   - Terms: show the links from `setup.py terms <skill>`, ask "did you read and accept them?", and only on
     yes run `install setup.licence:<skill>`. One acceptance per skill.
   - Key: `install setup.key:ALPHAGENOME_API_KEY` opens `~/.env` in the text editor with an empty line for
     the key. Tell them where to get one (https://deepmind.google.com/science/alphagenome/), to paste it
     after `=`, save, close and say "pronto"; then `setup.py key` confirms it is set.
7. **Verify**: `uv run SETUP/scripts/setup.py verify`. It runs each analysis skill's `check` and, when
   the skill's reference study is bundled, its offline `reproduce` (seconds once the packages are
   downloaded; the first run downloads them). Show "ready" or what is still missing.
8. **Report**: `uv run SETUP/scripts/setup.py report` writes `environment-report.md` in the workspace.
   Tell the researcher it is there to keep, or to attach when asking for help, and give one example
   request for each installed analysis skill (e.g. "quais variantes do gene X mudam a expressão no tecido
   Y e não em Z?").

Copy this checklist into the conversation and tick it as you go:

```
alphagenome-skills-setup progress:
- [ ] explained; uv present (installed with approval if needed)
- [ ] agent app asked; diagnose shown
- [ ] each missing item asked one at a time and installed only on approval
- [ ] plan run again after the skills; terms shown and accepted one by one
- [ ] key set by the researcher in ~/.env (never in the chat)
- [ ] verify passed (or what is missing explained)
- [ ] environment-report.md written and pointed to
```

## Commands

| Command | Does |
|---|---|
| `setup.py diagnose [--json] [--no-network]` | System and architecture, uv, Node.js, skills, terms, keys (presence only), disk, network |
| `setup.py plan [--json]` | What is missing, as questions with options, effects and the recommended one |
| `setup.py terms <skill>` | Links to the terms of use a dependency skill asks for |
| `setup.py install <key> --approved <key>` | `setup.uv`, `setup.node`, `setup.skill:<name>` (npx skills add, copy mode), `setup.licence:<skill>`, `setup.key:<VAR>` |
| `setup.py key [VAR]` | Whether a key is set; never its value |
| `setup.py verify [--json]` | `check` of each analysis skill, and its `reproduce` when the reference is bundled |
| `setup.py report` | `environment-report.md` in the workspace |

Global options (before the subcommand): `--workspace DIR`, `--target expression-analysis|splicing-analysis`
(repeatable; default: those installed, else expression-analysis), `--agents ID…`.

## Common Mistakes

- **Installing everything at once, or on a general "yes".** Each item gets its own question and approval.
- **Asking for the key in the chat, or reading `~/.env`.** The researcher pastes it in the editor; scripts
  check presence only.
- **Using `sudo` or administrator rights.** Not needed; on a server without them, see `references/platforms.md`.
- **Running `npx skills add` by hand without `--copy -y --agent …`.** It then asks questions the agent's
  shell cannot answer, or links instead of copying (links need developer mode on Windows).
- **Forgetting to run `plan` again after installing skills.** The terms of use come from the installed skills.
- **Overwriting an analysis skill.** `install setup.skill:<name>` does nothing when the skill is already there.

## Outputs

```
<workspace>/
  .claude/skills/, .agents/skills/   # skills installed with npx skills add (copy mode), per agent
  skills-lock.json                    # written by npx skills
  .licenses/                          # one notice per accepted terms of use
  .alphagenome-setup/state.json       # what this setup installed (source, hash, date) and the last verify
  environment-report.md               # the report for the researcher
~/.env                                # the keys, written by the researcher only
```

Read `references/platforms.md` for system-specific details and `references/troubleshooting.md` when a step fails.

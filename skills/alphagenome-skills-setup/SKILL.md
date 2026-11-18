---
name: alphagenome-skills-setup
description: "Prepares a researcher's computer (Windows, macOS or Linux) to use the AlphaGenome analysis skills expression-analysis and splicing-analysis: checks what is there, then, with the researcher's approval for each item, installs uv and Node.js, installs the analysis and dependency skills with npx skills add, explains and opens the terms of use and records those the researcher accepts, guides them through creating an AlphaGenome account and key and opens the project's .env (kept out of git) for them to paste it without it reaching the chat, verifies with each skill's check (and its offline reproduce when the reference is bundled) and writes an environment report. Use it whenever someone wants to install, set up or prepare the environment for expression-analysis or splicing-analysis, when one of their checks reports missing skills, terms or keys, or when a researcher without terminal experience asks how to start using these skills."
metadata:
  version: 0.2.0
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
- **The key never passes through the chat, and never reaches git.** The researcher pastes it into the
  project's `.env` (the only place the analysis skills read keys from) in a text editor; the setup adds
  `.env` to `.gitignore` first, because the AlphaGenome
  terms make the key personal (not to be shared, even with colleagues). Scripts only check that it is
  there. Do not read, print, copy or `cat` a `.env`.
- **No administrator rights, no Docker.** Everything installs for the user; Windows may still ask for
  permission when Node.js is installed with winget.
- **The agent's shell cannot type into prompts**, so every command here runs without questions.
- **A skill is installed only in the folder the agent app reads.** Each app lists the skills of its own
  folders and no others (Claude: `.claude/skills/`; Antigravity, Codex, Cursor, Gemini CLI, OpenCode:
  `.agents/skills/`; or the app's folder in the home). A copy anywhere else, such as a `.skills/` folder or
  another app's folder, is never loaded, and the agent then answers with the dependency skills alone,
  deciding windows and tissues without the researcher. `diagnose` shows where each skill is for each app.

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
3. **Which agent app.** `plan` asks it first when no app is recorded, and nothing else is offered before
   the answer. On it: `setup.py install setup.agent:<id> --approved setup.agent:<id>` (once per app; ids in
   the `agents` table of `assets/targets.json`). The app is recorded in `.alphagenome-setup/state.json`, and
   every later command prepares the folders of the recorded apps; `--agents id[,id]` overrides it for one call.
4. **Diagnose**: `uv run SETUP/scripts/setup.py diagnose`. Show a short table: the app and the folder it
   reads, what is there, what is missing, copies kept where the app does not read them, anything unusual
   (a service that does not answer, little disk, a `.env` inside the project, Windows on ARM). Nothing is
   installed here.
5. **Plan**: `uv run SETUP/scripts/setup.py plan` (or `--json`). For each item, one at a time: its
   title, the question, why it matters, what it is about, each option's label and what happens, and
   the recommended one. Follow its "For the agent" line first (e.g. show the terms with `setup.py terms`).
   Tell the researcher each "Notice" line (a service that cannot be reached now and what it affects): there
   is nothing to install for it, and it is usually the service, not their computer.
6. **Install what they approve**: `uv run SETUP/scripts/setup.py install <key> --approved <key>`, with
   the item's key; one call per approved item. A "not now" is fine: say what will not work without it
   and go on. A failure does not stop the rest: explain it in plain words (`references/troubleshooting.md`)
   and offer the alternative.
   - After Node.js, a new terminal session (or restarting the agent app) may be needed before `npx` works.
   - After the dependency skills are installed, run `plan` again: their terms of use only appear then.
   - A copy of an analysis skill, or of this setup skill, found outside the app's folder comes as
     `setup.place:<name>`: on yes it is copied into the app's folder (the old copy is left as it is). Never
     copy or move skills by hand; the setup records what it placed.
   - Terms: `setup.py terms <skill> --open` prints what the terms cover and opens the pages in their
     browser. Explain those points in plain words (who runs the service, what the terms allow and forbid,
     e.g. AlphaGenome: non-commercial research only, no clinical use, personal key), say the pages are what
     counts, give them time to read, ask "did you read and accept them?", and only on yes run
     `install setup.licence:<skill>`. One acceptance per skill.
   - Key: read `references/keys.md` and guide them one step at a time: a Google account, the key page
     (https://alphagenome.google/api), sign in, accept the AlphaGenome terms (explained as above), create
     and copy the key. `install setup.key:ALPHAGENOME_API_KEY` opens that page in the browser and the
     project's `.env` in the text editor with an empty line for the key, after adding `.env` to
     `.gitignore`. They paste the key after `=`, save, close and say "pronto"; `setup.py key` confirms it
     is set, never showing it.
7. **Verify** (required): `uv run SETUP/scripts/setup.py verify`. It checks that each analysis skill and
   this one are in the folder of every recorded app, runs each analysis skill's `check` and, when its
   reference study is bundled, its offline `reproduce` (seconds once the packages are downloaded), and lists
   the services that cannot be reached. Show "ready" or what is still missing. The setup is not finished
   before verify and the report: without them nobody knows whether the app will find the skill.
8. **Report**: `uv run SETUP/scripts/setup.py report` writes `environment-report.md` in the workspace.
   Tell the researcher it is there to keep, or to attach when asking for help, and give one example
   request for each installed analysis skill (e.g. "quais variantes do gene X mudam a expressão no tecido
   Y e não em Z?").

Copy this checklist into the conversation and tick it as you go:

```
alphagenome-skills-setup progress:
- [ ] explained; uv present (installed with approval if needed)
- [ ] agent app asked and recorded (setup.agent); diagnose shown
- [ ] each missing item asked one at a time and installed only on approval
- [ ] plan run again after the skills; terms explained, opened and accepted one by one
- [ ] account and key created with the researcher (references/keys.md); key pasted in the project's .env, .env ignored by git
- [ ] copies outside the app's folder placed (setup.place), or the researcher declined
- [ ] verify passed (or what is missing explained); unreachable services told
- [ ] environment-report.md written and pointed to
```

## Commands

| Command | Does |
|---|---|
| `setup.py diagnose [--json] [--no-network]` | System and architecture, uv, Node.js, the app(s) and their skill folders, skills (and copies outside those folders), terms, keys (presence only), disk, network |
| `setup.py plan [--json] [--no-network]` | What is missing, as questions with options, effects and the recommended one; notices for unreachable services |
| `setup.py terms <skill> [--open]` | What the terms of use cover, in plain words, with their links (opened in the browser with `--open`) |
| `setup.py install <key> --approved <key>` | `setup.agent:<id>` (records the app), `setup.uv`, `setup.node`, `setup.skill:<name>` (npx skills add, copy mode, for the recorded apps that lack it), `setup.place:<name>` (copies a copy found elsewhere into the apps' folders), `setup.licence:<skill>`, `setup.key:<VAR>` |
| `setup.py key [VAR]` | Whether a key is set in the project's `.env` (the only place used); never its value |
| `setup.py verify [--json] [--no-network]` | Each skill in the apps' folders, `check` of each analysis skill, its `reproduce` when the reference is bundled, unreachable services |
| `setup.py report` | `environment-report.md` in the workspace |

Global options (before the subcommand): `--workspace DIR`, `--target expression-analysis|splicing-analysis`
(repeatable; default: those installed, else expression-analysis), `--agents id[,id]` (default: the recorded
apps).


## Common Mistakes

- **Installing everything at once, or on a general "yes".** Each item gets its own question and approval.
- **Asking for the key in the chat, or reading a `.env`.** The researcher pastes it in the editor; scripts
  check presence only.
- **A `.env` that git does not ignore.** A commit would publish a personal key; `diagnose` and `key` warn,
  and `install setup.key:…` fixes `.gitignore` first.
- **Asking "do you accept?" before explaining.** Explain what the terms cover and open the pages first.
- **Using `sudo` or administrator rights.** Not needed; on a server without them, see `references/platforms.md`.
- **Running `npx skills add` by hand without `--copy -y --agent …`.** It then asks questions the agent's
  shell cannot answer, or links instead of copying (links need developer mode on Windows).
- **Forgetting to run `plan` again after installing skills.** The terms of use come from the installed skills.
- **Overwriting an analysis skill.** `install setup.skill:<name>` does nothing when the skill is already there.
- **Copying a skill into a folder by hand, or leaving it where it was downloaded.** The app only lists the
  skills of its own folder; `setup.place:<name>` copies it there and records it.
- **Ending the setup without `verify` and `report`.** Until verify passes, nobody knows whether the app finds
  the analysis skill.

## Outputs

```
<workspace>/
  .claude/skills/, .agents/skills/   # skills installed with npx skills add (copy mode), per agent
  skills-lock.json                    # written by npx skills
  .licenses/                          # one notice per accepted terms of use
  .env                                # the keys, pasted by the researcher only; listed in .gitignore
  .gitignore                          # gets a '.env' line
  .alphagenome-setup/state.json       # the agent app(s), what this setup installed or placed (source, hash, date) and the last verify
  environment-report.md               # the report for the researcher
```

Read `references/keys.md` before guiding the key, `references/platforms.md` for system-specific details and
`references/troubleshooting.md` when a step fails.

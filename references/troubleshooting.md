# When a step fails (alphagenome-skills-setup)

The message to look for, what it means, and what to say and do. Say it to the researcher in plain words.

| Message or symptom | Meaning | What to do |
|---|---|---|
| `uv: command not found` / `'uv' is not recognized` | uv missing, or installed but not on this session's PATH | Workflow step 2; after installing, use the full path or restart the agent app |
| `npx (Node.js) is not available` | Node.js missing, or not on PATH yet | `install setup.node`; restart the agent app after installing |
| `[SKILLS_INVALID_AGENTS] Invalid agents` | an `--agents` id npx does not know | use an id from `assets/targets.json` (`agents`) or the list npx prints |
| `<skill> was not installed (exit …)` with `Repository not found` / `Authentication failed` | the analysis skill's repository is private and the researcher is not signed in to GitHub with access | ask whether they have access; if so, `gh auth login --web` and `gh auth setup-git` (they sign in in the browser), then install again; if not, they ask the maintainers for access |
| npx hangs or `ETIMEDOUT` / `ECONNRESET` | network, proxy or firewall | `diagnose` (network lines); `references/platforms.md`, university network |
| `EPERM` / `operation not permitted` on Windows | antivirus or OneDrive locking a file | try again; move the workspace out of OneDrive; allow it in the antivirus |
| `setup.py key` says `not set` after the researcher saved | the key line was not saved, or was pasted with spaces, quotes or a line break | run `install setup.key:…` again (it reopens the project `.env`); the line must be `ALPHAGENOME_API_KEY=<key>` with nothing else; save |
| `verify`: `check` lists missing terms | a dependency skill's terms were not accepted in this workspace | `plan`, then the terms item |
| `verify`: `check` lists a missing skill | installed for another agent folder, or not installed | `diagnose` shows where each skill is; install it with the right `--agents` |
| `verify`: `unknown feature flag: 'sse3'` | an analysis skill without the pinned polars, on Windows on ARM | update the analysis skill to a version with pinned dependencies |
| `verify`: `reproduce` slow the first time | it downloads the Python packages (a few hundred MB) once | wait; later runs take seconds |
| `verify`: reproduce `skipped` | this copy of the analysis skill has no reference study | nothing to fix; `check` passing is enough |
| OLS or Ensembl calls that never finish (during an analysis) | a rate-limit reservation left far in the future by an interrupted run | the analysis skill's `fetch` resets its own; outside it, restart the computer or remove the old `polite-http` lock files in the temporary folder |

When nothing here fits: run `setup.py report`, read `environment-report.md` with the researcher, and
explain the failing line; do not try random fixes on their computer.

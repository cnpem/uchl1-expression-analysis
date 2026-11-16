# Platform notes (alphagenome-skills-setup)

What differs between systems, and what to tell the researcher. Read the section of their system when
`diagnose` shows something unusual or an install fails.

## Contents

- [Windows 10/11](#windows-1011)
- [macOS](#macos)
- [Linux](#linux)
- [University or company network](#university-or-company-network)

## Windows 10/11

- **uv**: `powershell -NoProfile -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"` installs
  it in `%USERPROFILE%\.local\bin` without administrator rights. The current session may not see it on
  PATH: use the full path, or restart the agent app.
- **Node.js**: with winget (present on current Windows 10/11), `setup.py install setup.node` runs
  `winget install -e --id OpenJS.NodeJS.LTS`; Windows may show a permission prompt, which the researcher
  accepts on screen. Without winget, the nodejs.org download page opens; the researcher runs the `.msi`.
  Restart the agent app afterwards so `npx` is on PATH.
- **Windows on ARM (Snapdragon and similar)**: `diagnose` shows `native: ARM64`. Python under emulation
  works; the analysis skills pin polars 1.31.0, the last version that runs there. Newer polars fail with
  `unknown feature flag: 'sse3'`.
- **Paths with spaces, accents or OneDrive** (`Área de Trabalho`, `OneDrive\Documentos`): they work, but
  synchronised folders can lock files while a study writes. Prefer a workspace outside OneDrive.
- **Antivirus** can block the uv installer or quarantine `uv.exe`: the researcher allows it in the
  antivirus, or IT does.
- **Console text**: the scripts write UTF-8; accents look right in the agent app.

## macOS

- **uv**: `curl -LsSf https://astral.sh/uv/install.sh | sh` installs it in `~/.local/bin`. New terminals
  see it; the agent may need the full path in the current session.
- **Node.js**: with Homebrew, `brew install node`; otherwise the nodejs.org page opens and the researcher
  runs the `.pkg` installer (it asks for the Mac password, which the researcher types in the installer,
  never in the chat).
- **Apple Silicon and Intel** both work.
- **Gatekeeper** may ask once to allow a downloaded program; the researcher allows it in System Settings.

## Linux

- **uv**: `curl -LsSf https://astral.sh/uv/install.sh | sh` (or `wget -qO- … | sh` without curl).
- **Node.js without administrator rights**: fnm installs it for the user:
  `curl -fsSL https://fnm.vercel.app/install | bash`, open a new shell, then `fnm install --lts`.
  With administrator rights, the distribution's package (`apt install nodejs npm`, `dnf install nodejs`)
  is fine, but this skill never runs `sudo`: the researcher (or IT) does it.
- **Servers without a desktop**: `install setup.key:…` cannot open an editor; the researcher opens
  `~/.env` with `nano ~/.env`, pastes the key after `=`, saves (Ctrl+O, Enter) and exits (Ctrl+X).
- **Old servers** with an old glibc may not run uv's Python; tell the researcher to ask IT, or to use a
  newer machine.

## University or company network

- `diagnose` tests the domains the setup and the analysis need (GitHub, npm, PyPI, Ensembl, EBI). An
  "unreachable" one usually means a proxy or firewall.
- **Proxy**: set `HTTPS_PROXY` (and `HTTP_PROXY`) for the session, and `npm config set proxy …` for npx; IT
  gives the address.
- **Institution certificates** (TLS inspection): `UV_NATIVE_TLS=1` makes uv use the system's
  certificates; for npm, `NODE_EXTRA_CA_CERTS` pointing to the institution's certificate file.
- Several studies at once on one computer share the Ensembl and OLS rate limits; the analysis skills
  handle it, but they run slower.

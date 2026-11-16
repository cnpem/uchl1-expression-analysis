# Getting the access keys (alphagenome-skills-setup)

How to guide a researcher who has never had an API key. Say it in their language, one step at a time,
and wait for "pronto" (or "done") before the next one. Steps checked on 2026-11-16; if a page looks
different, follow what it shows and tell the researcher so.

## AlphaGenome key (`ALPHAGENOME_API_KEY`, required)

What to say first: the analysis asks AlphaGenome (Google DeepMind) for predictions, and AlphaGenome
needs to know who is asking. The key is free for non-commercial research, **personal** (the terms forbid
sharing it, even with colleagues), and it goes into a file only on their computer, never into the chat.

1. **Google account.** They need one (a Gmail or an institutional Google Workspace account). Without one:
   https://accounts.google.com/signup.
2. **Open the key page**: `setup.py install setup.key:ALPHAGENOME_API_KEY --approved setup.key:ALPHAGENOME_API_KEY`
   opens https://alphagenome.google/api in the browser and the `.env` file in the text editor. (The page
   is also reachable from https://www.alphagenomedocs.com, *Get an API key*.)
3. **Sign in** on that page (*Sign in*, top right) with the Google account.
4. **Terms.** The page asks them to accept the AlphaGenome terms. Before they accept, explain the main
   points with `setup.py terms alphagenome-single-variant-analysis`: non-commercial use only, no clinical
   use, the key is personal, results can be published with the notice the Output Terms ask for.
5. **Create the key** and **copy** it (the page shows a button to copy it).
6. **Paste it in the `.env` file** that opened: right after `ALPHAGENOME_API_KEY=`, on the same line, with
   no spaces or quotes. Save (Ctrl+S, or Cmd+S on a Mac) and close.
7. They say "pronto". Run `setup.py key`: it answers `set` or `not set`, never the key.

If `not set`: the line was not saved, or has spaces, quotes or a line break; reopen the file with the same
install command and check the line. Do not ask them to read the key to you.

## NCBI key (`NCBI_API_KEY`, optional, splicing-analysis only)

It only raises the NCBI request limit (dbSNP, ClinVar); without it the splicing analysis is slower.

1. NCBI account: https://www.ncbi.nlm.nih.gov/account/ (they can sign in with Google or ORCID).
2. Account settings: https://www.ncbi.nlm.nih.gov/account/settings/, section *API Key Management*,
   *Create an API Key*, copy it.
3. `setup.py install setup.key:NCBI_API_KEY --approved setup.key:NCBI_API_KEY` and paste it after
   `NCBI_API_KEY=`, as above.

## Where the keys live

Only in the project's `.env` (the project folder): the analysis skills read keys from nowhere else (not
`~/.env`, not environment variables), so each study runs with its own project's key. The setup adds `.env`
to the project's `.gitignore` so it is never committed, since the AlphaGenome key is personal. A key left
in `~/.env` is not used; the setup says so. Never read, print, copy or attach a `.env`.

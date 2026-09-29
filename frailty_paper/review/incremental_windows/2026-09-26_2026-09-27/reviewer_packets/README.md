# Supplemental PubMed title/abstract screening

This packet contains seven unique records retrieved for publication dates 2026-09-26 through 2026-09-27. It is a supplemental queue and is not merged into the frozen 56,117-record canonical queue.

Two people must screen independently. Each reviewer should receive only their own ZIP and must not compare decisions before both returns are collected. Read the included frozen `protocol.md` before screening. In the single TSV batch, preserve `screening_id` and every source-information column; fill only `decision`, `reason`, `reviewer`, and optional `notes`. Use `include`, `exclude`, or `unclear` for `decision`; identify yourself consistently in `reviewer`. Do not make full-text decisions in this file.

Return the completed TSV without renaming it, together with this README and `protocol.md`. Keep a copy of the original blank archive. The project team will validate both returned reviewer packets before applying decisions to this supplemental queue. Do not merge these records into the canonical queue automatically.

## Return validation and import

From the repository root, copy the original blank packet directory to a staging directory, then extract each return into its matching reviewer slot. Keep the original directory and ZIPs unchanged:

```sh
cp -R frailty_paper/review/incremental_windows/2026-09-26_2026-09-27/reviewer_packets \
  frailty_paper/review/incremental_windows/2026-09-26_2026-09-27/reviewer_packets_returned
unzip -o reviewer_1_return.zip -d \
  frailty_paper/review/incremental_windows/2026-09-26_2026-09-27/reviewer_packets_returned
unzip -o reviewer_2_return.zip -d \
  frailty_paper/review/incremental_windows/2026-09-26_2026-09-27/reviewer_packets_returned
```

Dry-run each reviewer independently first. Both packets must pass before applying either slot:

```sh
.venv/bin/python frailty_paper/scripts/56_import_pubmed_window_screening.py \
  --packet-dir frailty_paper/review/incremental_windows/2026-09-26_2026-09-27/reviewer_packets_returned \
  --slot reviewer_1
.venv/bin/python frailty_paper/scripts/56_import_pubmed_window_screening.py \
  --packet-dir frailty_paper/review/incremental_windows/2026-09-26_2026-09-27/reviewer_packets_returned \
  --slot reviewer_2
```

Then repeat each command with `--apply`. The importer verifies the frozen protocol and packet manifest hashes, full record coverage, immutable source rows, decision values and distinct reviewer identities. It refuses incomplete or conflicting returns, creates a hash-named byte-preserving queue backup, and never edits the canonical 56,117-record queue. It reports disagreements but does not adjudicate them.

The verified blank supplemental archives are `reviewer_1_2026-09-26_2026-09-27.zip` (SHA-256 `3762dd1b930085e5b37bb57c9e31fbeae1401589c48820e6989a576f3076aaa7`) and `reviewer_2_2026-09-26_2026-09-27.zip` (SHA-256 `00ef5e17a2768912ef5ef9090bcaf1f8c7bc39ec5c609e66fe4813c2af1cdac6`). Each contains one seven-record batch for only its assigned reviewer; both ZIP integrity checks passed. The importer also pins the blank queue SHA-256 (`7461d548b509d5a035c2798e211d9b90d464018861a9193be4f23b56bbe03ba7`), frozen protocol SHA-256 (`1013526fa40de6708ed9e547e09c1cb0ae3de0c27d8cbe2d5d04a99ffa796abb`) and packet manifest SHA-256 (`83989c4c13b4354940dd2b125a0afbfde886c89659b0c393cb4604efd859f8a4`). It uses the original queue hash-named backup to verify baseline provenance after the first slot has been imported.

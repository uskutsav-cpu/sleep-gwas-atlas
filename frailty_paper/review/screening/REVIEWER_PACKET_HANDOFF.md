# Independent title/abstract screening packet handoff

## What is ready

The frozen title/abstract queue contains 56,117 records and has zero decisions. Two independent packet sets are prepared under `reviewer_packets_20260926T0835Z/`: `reviewer_1/` and `reviewer_2/`. Each set contains 113 TSV batches (112 batches of 500 records and one final batch of 117). The manifest records the source queue and protocol hashes and per-batch content hashes. Packet generation did not edit the canonical queue or make screening decisions.

## Shareable reviewer archives

Separate compressed archives are prepared for handoff. Each archive contains only that reviewer’s 113 TSV batches, a reviewer-specific manifest and instructions; do not send both archives to the same reviewer.

- `~/Documents/Codex/frailty-review-handoff-2026-09-26/reviewer_1_title_abstract_screening.zip` — SHA-256 `88bfc194c6c900e87f335d67a2dbda2c3e04daf65f6110e7a009ec9c1ddf9d41`
- `~/Documents/Codex/frailty-review-handoff-2026-09-26/reviewer_2_title_abstract_screening.zip` — SHA-256 `0b87aaebfdea35022ccf59aefb648d06f40e5ece665b2e9921fa35a19b42f4f6`

The archives are outside the Git repository. The outbound archive SHA-256 values and instructions are also recorded in `~/Documents/Codex/frailty-review-handoff-2026-09-26/README.md`.
Both archives were extracted into a temporary staging copy and passed the full importer dry-run for their respective slots (56,117 records each; zero decisions applied). Evidence: `analysis/reviewer_archive_validation_2026-09-26_2102.md`.

The user has confirmed that two distinct reviewers are ready. Assign Reviewer 1 the reviewer 1 archive and Reviewer 2 the reviewer 2 archive; never send both full-size reviewer archives to one person. There is also a separate seven-record PubMed supplement for 2026-09-26–27. Give each person only their corresponding supplemental ZIP, in addition to their canonical queue archive:

- `review/incremental_windows/2026-09-26_2026-09-27/reviewer_packets/reviewer_1_2026-09-26_2026-09-27.zip` — SHA-256 `3762dd1b930085e5b37bb57c9e31fbeae1401589c48820e6989a576f3076aaa7`.
- `review/incremental_windows/2026-09-26_2026-09-27/reviewer_packets/reviewer_2_2026-09-26_2026-09-27.zip` — SHA-256 `00ef5e17a2768912ef5ef9090bcaf1f8c7bc39ec5c609e66fe4813c2af1cdac6`.

Each supplemental ZIP contains only its assigned reviewer's seven-record batch plus protocol/instructions. Both passed ZIP integrity checks. Their return staging, dry-run and import procedure is in `review/incremental_windows/2026-09-26_2026-09-27/reviewer_packets/README.md`; use `scripts/56_import_pubmed_window_screening.py`. This supplemental queue remains separate from the frozen 56,117-record corpus and its PRISMA counts.

## Reviewer instructions

Assign one complete directory to each distinct human reviewer. Reviewers must work independently and should not exchange decisions before both packet sets are returned. Keep `screening_id` and every source-information column unchanged. In each row, complete only:

- `decision`: `include`, `exclude`, or `unclear`;
- `reason`: concise protocol-based rationale when useful;
- `reviewer`: the reviewer's name or agreed unique identifier;
- `notes`: optional note for context.

Do not enter full-text decisions in these files. Records needing full text proceed to the separate full-text stage after title/abstract screening and reconciliation.

## Return and import

Return each reviewer’s complete 113-file directory without renaming the batch files. If a reviewer returns a ZIP, preserve it and extract it into a separate staging copy; never extract over the original blank packets. For example, from the repository root:

```sh
cp -R frailty_paper/review/screening/reviewer_packets_20260926T0835Z \
  frailty_paper/review/screening/reviewer_packets_20260926T0835Z_returned
unzip -o ~/Documents/Codex/frailty-review-handoff-2026-09-26/reviewer_1_title_abstract_screening.zip \
  -d frailty_paper/review/screening/reviewer_packets_20260926T0835Z_returned
```

For reviewer 2, extract their archive into the matching reviewer slot of a staging copy too. You may use one staging copy containing both returned slots or separate copies per reviewer; keep the original blank packet tree intact. Keep the copied root `packet_manifest.json` and both slot directories, replacing only the matching reviewer slot’s batches with returned files. Then validate each slot independently without changing the canonical queue:

```sh
.venv/bin/python frailty_paper/scripts/import_independent_screening_packets.py \
  --packet-dir frailty_paper/review/screening/reviewer_packets_20260926T0835Z_returned \
  --slot reviewer_1

.venv/bin/python frailty_paper/scripts/import_independent_screening_packets.py \
  --packet-dir frailty_paper/review/screening/reviewer_packets_20260926T0835Z_returned \
  --slot reviewer_2
```

Only after both dry-run imports validate should the reviewed slot be merged by rerunning the corresponding command with `--apply`. The importer checks complete queue coverage, immutable source rows, valid decisions, reviewer identity, and cross-reviewer independence, and creates a byte-preserving queue backup before its atomic update. Then validate the combined queue:

```sh
make -C frailty_paper validate-review
```

The 56,117 records in this queue are not a substitute for licensed database searches. Embase, Scopus, Web of Science, and PsycINFO exports still need to be obtained and imported under the documented source-preservation workflow. No packet decisions or licensed exports are present yet.

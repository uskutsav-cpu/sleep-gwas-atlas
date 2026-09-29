# Brain6 analysis artifacts

This directory contains the checked global-map extraction and its current evidence/coverage package. It does not replace or mutate the atlas or the protected Track B analyses.

The current integrated scientific snapshot and evidence-linked readiness assessment is [`FINAL_BRAIN6_REPORT.md`](FINAL_BRAIN6_REPORT.md); its report provenance is recorded in `FINAL_BRAIN6_REPORT.provenance.json`. The report is explicitly provisional while the separate roundoff LAVA family continues.

## Reproducible entry point

Run `make -C brain6 help` from the repository root to see supported targets. `make -C brain6 check` sequentially rebuilds the interim sensitivity ledger and software Table S17, runs both test suites, validates current outputs, and audits available production receipts. `make -C brain6 audit-checkpoints` runs the read-only receipt check alone; it verifies every available PLACO chunk and LAVA locus receipt/output and reports whether the run-manifest counters match. A mismatch can indicate a stale run manifest; check the relevant run's final decision before taking action. Global rebuilds run a fast `analysis-deps` preflight and use a writable temporary Matplotlib cache. Analysis targets remain explicit; there is no umbrella target that starts long production jobs without the user selecting them. Override mounted-volume locations with `LAVA_INPUT_ROOT`, `LAVA_OUTPUT_ROOT`, and `RUNTIME_RSCRIPT` when the external research volume is mounted at another path.

## Rebuild the global map

From the repository root, install the analysis extra into the Brain6 virtual environment, then use the default Makefile interpreter:

```bash
extensions/brain6/.venv/bin/python -m pip install -e 'extensions/brain6[analysis,test]'
make -C brain6 global
```

The builder reads `config/analysis_panel.tsv`, `results/atlas/traits.tsv`, `results/atlas/trait_pairs.tsv`, and `results/analysis/phase1_master_analysis.tsv`. It requires the exact 12×6 pair family, validates original 396-family q values, refuses to overwrite outputs, and computes no new significance correction. It writes the 72-row TSV, its figure-source table, the family config/checksum, five figures in PNG/PDF, and descriptive profile similarity.

`python brain6/scripts/build_profile_spearman.py` adds a rank-based descriptive sensitivity for the six disorder profiles. It reads the locked family dimensions and global map, writes a separate non-overwriting table plus provenance, and makes no inferential claims or new discovery family.

The `analysis` extra adds pandas for this builder; NumPy, SciPy, and Matplotlib are already package runtime dependencies. Output scripts use a separate output namespace under `brain6/`.

Use `python brain6/scripts/build_locked_global.py --refresh-figures` to re-render only the figures after the script confirms the existing 72-row table, similarity table, and family config still match the locked atlas. `python brain6/scripts/build_headline_claim_sensitivity.py` rebuilds the interim, source-hashed sensitivity ledger from current global, replication, PLACO, candidate-locus, and LAVA evidence. Untested sensitivities remain labeled; the interim ledger is not a final robustness assessment. `python brain6/scripts/validate_current_outputs.py` checks the table values against the original atlas, follow-up SHA, figure inventory, supplementary tables, candidate-locus artifacts, and interim sensitivity evidence hashes.

`make -C brain6 lava-not-run-causes` builds a checksum-bound row for each of the 3,720 frozen canonical v3 `NOT_RUN` cells. `make -C brain6 lava-not-run-summary` adds immutable seven-trait and 2,495-locus breakdowns for Tables S25–S26, bound to the complete aggregate, decision, latest receipt audit, and cell-level causes. `make -C brain6 cross-layer-evidence` builds the five-primary-pair evidence integration from the inherited global map, replication ledger, partial-family PLACO outputs, and frozen LAVA decision. These outputs are diagnostic/descriptive and do not change the failed-QC decision. `make -C brain6 validate` rebuilds or verifies them before validating current outputs.

## Frozen deep-follow-up selection

`scripts/freeze_followup_manifest.py` reproduces the already-reviewed post-atlas pair lock from `extensions/brain6/work/overnight-v03/`. It refuses to write over an existing v1 table and stores a SHA256 beside `config/deep_tracks_v1.tsv`. It records the five existing selected pairs and Alzheimer's `NO_ELIGIBLE_PAIR`, preserving unknown source and overlap facts.

Run the Brain6 extension tests from `extensions/brain6/` so its local `brain6` package is imported rather than the root-level results directory:

```bash
cd extensions/brain6
PATH="/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin:$PATH" \
  .venv/bin/python -m pytest -q -p no:cacheprovider tests
cd ../..
python -m pytest -q brain6/scripts/tests
python brain6/scripts/validate_current_outputs.py
```

Put the checksum-validated R 4.3.3 runtime first on `PATH` so native adapter tests probe the same runtime used by production. The current complete extension suite passes 308 tests; three native checks are skipped because `susieR`, `TwoSampleMR`, and PLINK 1.9 are unavailable. The separate Brain6 script suite passes 111 tests. The output validator checks the current published artifact hashes and explicitly reports that no final freeze exists.

This is a post-atlas selection. The selection file is not prospective and does not correct the original global screen's selective-inference properties.

## Evidence boundary

`STATUS.md`, `BLOCKERS.md`, `paper/REVIEWER_2_AUDIT.md`, and `FINAL_READINESS_REPORT.md` record which claims are supported, unavailable, or not justified. Do not treat readiness files, synthetic tests, or plans as empirical variant-level results. The final release freeze is deferred until the required real-data analyses and validations are complete.

## Dense-input and local-rg provenance

The seven selected dense GWAS inputs and 21 preparation receipts are verified; provenance and source limitations are recorded in `manifests/locked_dense_input_audit.tsv`, `qc/source_archive_audit.tsv`, and the readiness report. The sealed UK Biobank v1.1 reference and pinned LAVA 0.1.5 runtime were used for the corrected canonical family. Source-card hashes, participant-overlap limitations, and nonuniform INFO coverage remain part of the scientific interpretation.

### Current LAVA state

The pair-independent canonical v3 run is complete under immutable run ID `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`, with its output in `work/lava-canonical-v3-production/`. Its frozen execution policy is four workers and five loci per persistent chromosome batch. The verifier checked all 2,495 locus receipts and 17,465 canonical trait-by-locus cells: 13,745 `TESTED`, 3,720 `NOT_RUN`, zero `FAILED`, and zero invalid receipts. Since the frozen maximum is 873 untested cells (5%), the immutable family decision is `FAILED_QC_NOT_PROMOTED`. Do not promote local-rg results or launch its bivariate stage.

The earlier v2 family remains preserved as a failed aggregation. The separate roundoff run completed all 2,495/2,495 locus receipts, but its aggregation failed on a repeated insomnia/locus-2 univariate disagreement; it produces no family decision or promoted local-rg result. A read-only audit found zero receipt errors, and the 12,475-slot comparison plus the isolated locus-740 retry are recorded separately. The retry returned `UNIVARIATE_UNDERPOWERED` and was not merged into the run. Protected Track B remains unavailable: the archived candidate lacks aligned-input/runner provenance and fails its terminal monitor gate, so it was not imported. See `qc/protected_track_b_archive_audit.md`, the roundoff comparison under `results/lava/`, and `results/project_state_inventory_20260924.json`. Never merge or copy receipts between runs. Read current state with:

```bash
make -C brain6 audit-checkpoints
python brain6/scripts/validate_current_outputs.py
```

The first command audits available receipts; the second verifies the frozen canonical v3 decision and the repository's current output bindings. The historical `run_lava_family.py` and `lava-roundoff-run` commands below describe earlier run variants and are not a route to repair or promote v2/v3 results. Any scientifically material new local analysis requires a separately reviewed design and a distinct run identity/output root.

The environmental retry has already been completed as `envfix-20260924-001`; its receipt-bound result is `UNIVARIATE_UNDERPOWERED`, and it remains isolated from the original roundoff family. Once complete receipt sets and the retry artifact exist, rebuild the descriptive comparison and Table S24 with `make -C brain6 lava-roundoff-compare`. The comparison target refuses to run unless both source families pass the full receipt audit.

The benchmark, receipt-bound runtime profile, aggregation-failure audit, final decision, and rationale are documented in `qc/lava_canonical_v3_benchmark_20260923.md`, `qc/lava_canonical_v3_runtime_profile_20260924T0731Z.md`, `qc/lava_family_aggregation_failure_audit_20260923.md`, and `BLOCKERS.md`.

## Publish completed PLACO+ v3 pair outputs

After a pair has all frozen chunks and a QC-passing collation, publish its full variant table and provenance with:

```bash
python brain6/scripts/publish_placo_results.py --pair insomnia__mdd
```

The publisher verifies every chunk receipt, the collation receipt, row denominator, SNP uniqueness, P/q ranges, numerical-failure count, and the frozen headline threshold. It writes `results/placo/<pair>/variants.tsv.gz`, the source receipt and QC status, and updates `results/placo/placo_master.tsv`. Publishing is idempotent for identical bytes and refuses to replace a different existing result. Its focused tests are in `brain6/scripts/tests/test_publish_placo_results.py`.

## Partial PLACO candidate-locus grouping

`config/shared_locus_rule_v1.json` freezes the preannotation candidate rule: pair-QC-passed PLACO variants below the five-track threshold are greedily clumped against LAVA's UKB v1.1 LD with `r² ≥ 0.1` within 500 kb. The extractor uses the pinned LAVA 0.1.5 `.bcor` reader and requires exact SNP ID, chromosome, and position matches to the reference. Run `python brain6/scripts/clump_placo_candidates.py` after more pair tables are published to refresh the explicitly partial candidate files under `results/loci/`. They are PLACO-only candidates, not final independent loci or shared biological mechanisms. `validate_current_outputs.py` checks the rule lock, result hashes, variant/locus accounting, and partial-family status.

Before any functional annotation, `config/shared_locus_evidence_tiers_v1.json` freezes the pair-by-locus evidence-tier criteria and binds them to the locked global map, five-pair LAVA family, five-track PLACO family, shared-locus rule, deep-track selection, and replication baseline by SHA256. It requires a complete QC-passing PLACO family and an eligible LAVA family decision; partial-family candidate regions cannot receive final tiers. No annotations were inspected to define the policy. No tiers are assigned because the five-track PLACO family is incomplete and canonical LAVA v3 failed its frozen QC gate.

`python brain6/scripts/validate_current_outputs.py` also validates the tier-policy checksum and every source binding, and fails if a final `shared_loci.tsv` appears before those prerequisite families are complete.

The four available pair-QC-passed PLACO candidate tables are visualized in `results/placo/figures/figS_partial_placo_candidate_intervals_v1.pdf` (300-dpi PNG preview alongside it). The supplementary figure is explicitly partial because protected insomnia–ADHD Track B is unavailable; its intervals are PLACO-only candidates, not final independent loci or causal mechanisms. `make -C brain6 placo-candidate-figure` verifies the existing checksum-bound outputs and refuses to overwrite them.

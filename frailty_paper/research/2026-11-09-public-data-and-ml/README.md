# Sleep–Frailty public-data and molecular-model research snapshot

Snapshot: **9 November 2026**. This directory preserves this research workspace's scripts, reports, tables, manuscript materials, and current native-model evidence. It supplements the existing repository analyses without replacing their results.

**The project has not produced a validated clinical prediction model. The completed native-model comparison failed its prespecified accuracy gate. The separate permitted AlphaGenome-related input experiment has not started scientific fitting.**

## Current scientific findings

| Evidence | Verified status | Interpretation |
| --- | --- | --- |
| Native DNA-model comparison | Eight full fits; all 24 saved effect-magnitude results independently checked; matched initialization, examples and schedules verified | Updating the DNA reader did not meet the required accuracy improvement across the three omitted neural studies. |
| Common-pass decision | Pass 2, selected by the original three-source rule | Mean adapted skill was −0.0005948435: average relative error was approximately 0.0595% worse than the typical-effect baseline. The gain over the matched frozen control was about 0.11285 percentage points, below the required 0.5 points. |
| Nearby and cell-specific differences | Saved pair results checked at the selected pass | Basic preservation checks passed, but these results do not establish a useful predictor across experiments. The Retallick source contains only nine nearby pairs and no cross-cell pairs. |
| Confidence prototype | The earlier archived-predictor experiment qualified in 0/7 source checks | No useful confidence policy was released. This is a result for that historical prototype, not a new confidence evaluation of the native model. |
| Permitted composite AVI inputs | Allele-query corrections, train-only scaling and numerical input/model checks preserved | No predictive benefit has been demonstrated. Real execution, continuation and resource checks of the separate trainer remain pending. Any composite-score benefit cannot be assigned uniquely to AlphaGenome. |
| Untouched-study test | Its outcomes remain excluded from model selection | No new-study benefit, clinical validation or independent genetic-score validation is claimed. |

The native comparison required both at least 1% mean magnitude skill over the typical-effect baseline and at least a 0.5-percentage-point gain over the frozen control. Both requirements failed. Its other four whole-source roles were therefore not run. No learning-rate, architecture or epoch extension rescues this finite comparison.

The model predicts measured molecular effects in reporter experiments. It does **not** predict a person's insomnia or frailty risk. Historical AlphaGenome/API reports and exploratory model versions are preserved as historical evidence; ordinary API/Atlas predictions are not teacher labels or added inputs in the new permitted-score experiment.

## Start with these records

- [Independent replay of the native comparison decision](outputs/sleep_frailty_ml_resumable/native_comparison_independent_decision_audit.json)
- [Original fixed comparison plan](outputs/sleep_frailty_ml_resumable/finite_stable_native_plan.json)
- [Completed Retallick matched comparison](outputs/sleep_frailty_ml_resumable/completed_Retallick2026_matched_comparison.json)
- [Completed Kosicki matched comparison](outputs/sleep_frailty_ml_resumable/completed_Kosicki2025_matched_comparison.json)
- [Completed Salomon matched comparison](outputs/sleep_frailty_ml_resumable/completed_Salomon_matched_comparison.json)
- [Confidence prototype decision](outputs/sleep_frailty_ml_confidence_revision/confidence_method_decision.json)
- [Untouched-study input-benefit protocol](outputs/sleep_frailty_ml_resumable/AVI_preparation/fresh_study_benefit_protocol.json)
- [Current plain-language model status](outputs/sleep_frailty_ml_resumable/PLAIN_LANGUAGE_PROGRESS.txt)
- [NHANES manuscript audit and readiness assessment](outputs/nhanes_manuscript_audit/audit_and_readiness.md)
- [Manuscript incorporating the observational NHANES sensitivity](outputs/nhanes_manuscript_audit/manuscript_with_nhanes.md)

## Keep observational and genetic evidence separate

The NHANES sensitivity concerns an **adapted deficit-accumulation Frailty Index** after removing direct sleep and fatigue items. It is distinct from the original Atkins FI and from the Fried phenotype. Its observational association does not establish genetic correlation, causation or independent genetic-score validation. The original genetic analyses and their limitations remain preserved separately.

The remaining scientific gaps include a genuine matched sleep/exhaustion-excluded frailty GWAS, independent genetic replication, stronger causal-variant-to-gene evidence, relevant molecular/cellular validation, and useful prediction confidence across experiments. This snapshot is not a submission-readiness claim.

## What is versioned

`work/` contains preserved analysis scripts and selected configuration/provenance files. `outputs/` contains reports, manuscript materials, derived tables, numerical audits, and the current native fit's checkpoints, predictions and recovery evidence. [snapshot_manifest.tsv](snapshot_manifest.tsv) binds copied files to their original bytes. [PUBLICATION_RECEIPT.json](PUBLICATION_RECEIPT.json) records the publication scope.

Large downloads, older binary model artifacts, and re-fetchable response caches remain external. Their retained paths, sizes and SHA-256 hashes are listed in [external_artifacts.tsv](external_artifacts.tsv). Original source/download receipts document reacquisition where available. This is an inventory, not a claim that every external file can be downloaded automatically. Installed environments, credentials, and local-only account/destination planning records are not distributed.

Preserved scripts may contain original machine paths and historical, unsuccessful methods. They are not all fresh-clone entry points. Do not interpret an archived engineering check or synthetic fixture as biological validation. The small replay below does not load the large DNA banks, train a model, use a GPU, or open reserved-study outcomes.

## Verify the release and replay the measured-effect results

```sh
python -m pip install -r requirements-snapshot.txt
PYTHONDONTWRITEBYTECODE=1 python verify_snapshot.py --replay-native
```

The verifier checks distributed file hashes and current-core Python syntax, then replays all 24 raw measured-magnitude metrics and the locked stopping decision in a temporary copy. It preserves the distributed evidence files. Passing these checks establishes snapshot integrity and numerical replay, not useful prediction accuracy or clinical readiness.

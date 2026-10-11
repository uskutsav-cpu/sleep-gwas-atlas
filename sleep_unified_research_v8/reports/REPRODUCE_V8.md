# Reproduce the V8 interpretation package

This package is additive to preserved V4/V7 results. It does not repeat the completed native campaign. Use Python 3 standard library and base R; no new software environment or controller is required.

From a package directory containing resource_tools, tables and source_review:

```sh
python3 resource_tools/build_interpretive_release.py --historical-release /path/to/preserved/V4/release_candidate --v8-root . --out release_candidate
python3 resource_tools/query_atlas.py --release release_candidate estimate --family core --sleep insomnia --outcome ibd
python3 resource_tools/query_atlas.py --release release_candidate estimate --family validation --sleep insomnia --outcome gwas_catalog_GCST90479148
Rscript resource_tools/plot_transport.R tables/insomnia_abdominal_pain_source_defined_estimates.tsv figures/insomnia_abdominal_pain_source_transport
Rscript resource_tools/plot_core_interpretive.R release_candidate/core_rg.tsv figures/core_interpretive_quarantine
```

The builder verifies the exact 36 held keys and source identities, refuses output overlapping the historical candidate, preserves all original inputs, withholds signed fields in a separate core view, and includes the reviewed exploratory result in a separate one-row family. The original seven other TSVs remain byte-identical. The manifest records candidate-input/output digests. Meaningful checks compare held/unheld rows and original testing fields rather than claiming that a manifest proves biological validity.

The historical V4 candidate can be rebuilt using its committed `sleep_unified_research_v4/cell_genomics_readiness_v1/scripts/build_candidate_release.py`; use a separate output location so no old candidate is overwritten. Required historical builder inputs and provenance are in the committed V4 handoff. Full native reruns require permitted SSD source access and the frozen environment/reference contract; this reporting rebuild does not.

The new native fit's exact plan, reused monitor/worker implementation hashes, frozen protocol, stock capture, logs, full-precision result and actual deletion vectors are preserved privately on the SSD. `../source_review/statistical_scientific_actual_result_review_v8.json` records independent arithmetic and input/result locations. No new executable controller was written. Do not rerun a completed fit merely to recreate this interpretation view.

The public code checkpoint excludes the candidate data directory, private human-scoring key, third-party source/paper extracts and all raw summary statistics. The complete candidate package is preserved locally under the output/SSD handoff, pending source-rights and human approval. Public script availability does not constitute cleared resource deposition.

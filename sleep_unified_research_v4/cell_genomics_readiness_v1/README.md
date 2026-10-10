# Sleep GWAS Atlas: resource-readiness package

This package prepares the completed Sleep GWAS Atlas V4 for an integrated resource manuscript. It addresses the seven priorities in the current goal while keeping the original 396-test core family, 1,200-test extension family, and 217-candidate external-validation family distinct.

The candidate machine-readable release is in `release_candidate/`. It contains allowlisted summary results, source-defined sleep-construct metadata and an outcome phenotype crosswalk for deterministic joins. It excludes GWAS bodies, local paths, source-body hashes and intermediate logs. It is a review candidate, not a public deposition: investigators must confirm source terms, metadata, permissions and final data-availability language before depositing it.

No manuscript prose is drafted here. The package contains evidence, definitions, use cases, figure-source mapping, reproducibility instructions and the remaining authoring decisions.

## Seven priorities

1. `01_resource_benchmark.md` and `benchmark_matrix.tsv` compare the nearest existing resources and bound the atlas contribution.
2. `02_source_harmonization.md` and `02_source_harmonization.tsv` document source-specific phenotype and processing metadata; unverified fields remain explicit.
3. `03_uncertainty_reliability.md` and `03_uncertainty_reliability_field_dictionary.tsv` distinguish reported marginal uncertainty from uncalibrated cross-estimate contrasts and unidentified measurement reliability.
4. `04_biological_use_cases.md` gives result-level examples of defensible source-aware queries without making a novelty or mechanism claim.
5. `release_candidate/` is the machine-readable summary-result candidate, schema and checksummed manifest.
6. `06_independent_usability_review.md` records the independent consumer review. The final candidate passed its defined packaging/usability scope; source-specific rights and metadata clearance still gate public deposition.
7. `07_reproducibility_and_authoring_handoff.md` gives the build/check commands, frozen-input inventory, figure map and critical path to human manuscript authoring.

`05_resource_release_review.md` explains the release scope and rights gate. `figure_selection.tsv` maps existing figures to result tables and interpretation limits. `08_execution_status_and_decision.md` records the exact critical path, optional provenance work, running-process qualifications and publication boundary.

## Build and validate

From the repository root, run:

```bash
python3 sleep_unified_research_v4/cell_genomics_readiness_v1/scripts/build_candidate_release.py
python3 sleep_unified_research_v4/cell_genomics_readiness_v1/tests/test_candidate_release.py
```

These commands read the already completed tables. They do not acquire sources, rerun LDSC, or modify frozen outputs. Python's standard library is sufficient.

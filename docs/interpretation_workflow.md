# Locked interpretation, graph, and robustness workflow

This workflow covers the remaining atlas interpretation layers without turning
missing data or a null analysis into a positive biological claim. Its scope is
frozen in `config/interpretation_analysis_policy.json`; external resources and
runtime families are enumerated in `config/interpretation_source_registry.tsv`.

## Pre-result source and task locks

`scripts/74_interpretation_preflight.py` verifies the exact 45-trait panel,
aligns the interpretation and downstream policies, audits every source, and
reports upstream readiness. A downloadable or licensed source is not ready
until its release, local path, exact byte count, and SHA256 have been curated.
Derived regulatory and QTL evidence is not ready until its canonical upstream
atlas artifact exists.

Only a production-ready preflight can reach
`scripts/75_prepare_interpretation_tasks.py`. That script freezes, before any
result is opened:

- every primary locus × six regulatory layers × four context domains;
- every one of 45 traits × five cell-type strategies × four domains;
- the four pathway resource families applied to high-confidence convergent
  genes only; and
- all 396 sleep/non-sleep pairs × two directions × five MR estimator families.

The task manifest records a deterministic input-scope hash and output path for
each unit. Existing task results cause locking to fail.

## Result recording and canonical evidence

Each external or derived analysis is imported through an explicit curator JSON
record and `scripts/76_record_interpretation_task.py`. A task ends as
`COMPLETED`, `NO_EVIDENCE_FOUND`, `ACCESS_BLOCKED`, or `NOT_APPLICABLE` with a
reason. `COMPLETED` requires a real exact-schema table; the other states cannot
smuggle in rows. Output and provenance are immutable rather than overwritten.

`scripts/77_collate_interpretation.py` revalidates every task and creates:

- `results/atlas/regulatory_elements.tsv`;
- `results/atlas/cell_types.tsv`;
- `results/atlas/pathways.tsv`;
- `results/atlas/causal_tests.tsv`; and
- `results/tables/interpretation_coverage.tsv`.

Cell-type and pathway P values receive BH correction within their frozen test
families. Only FDR-supported rows enter those two canonical evidence tables.
The complete null family remains auditable through coverage. Regulatory rows
must link a fine-mapped same-locus variant to an already supported gene; a
nearest-gene-only assignment is forbidden. MR rows need at least three
instruments, minimum F ≥ 10, and all required diagnostics. Robust causal
wording additionally needs corrected multi-estimator agreement, clear
pleiotropy diagnostics, correct Steiger direction, characterized or adjusted
sample overlap, and compatible signal-level colocalization.

Access-blocked tasks can be recorded for audit, but they cannot pass the final
interpretation or integrated-atlas gates.

## Atlas graph and robustness

`scripts/78_build_atlas_edges.py` creates deterministic, cross-linked graph
edges and a major-conclusion manifest. Only prespecified high-evidence classes
become major conclusions: primary global genetic correlations, cross-method
shared loci, convergent genes, multi-strategy cell types, FDR-supported
pathways, and robust causal-inference results.

Before sensitivity results are viewed,
`scripts/79_prepare_robustness_tasks.py` expands every major conclusion across
the nine locked robustness families. Applicability is frozen by conclusion
type. Methodologically non-applicable rows receive their prespecified reason;
applicable rows require a real checksum-bound sensitivity artifact.
`scripts/81_collate_robustness.py` and the independent
`scripts/53_validate_robustness.py` require the exact conclusion × family
matrix and reject unresolved result-changing contradictions.

## Current production blockers

The code and task contracts are ready, but production is not. At the current
repository state, only the already pinned GENCODE promoter resource is ready;
the other external releases still need exact pre-result curation, and four of
six upstream canonical inputs do not yet exist. No downstream task result has
been generated or claimed.

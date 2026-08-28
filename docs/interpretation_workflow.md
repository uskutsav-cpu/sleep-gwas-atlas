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
record and `scripts/76_record_interpretation_task.py`, except for the locked
GENCODE, SCREEN, HOCOMOCO, ABC, and promoter-capture Hi-C adapters, which produce the same
checksum-bound normalized record automatically. A task ends as
`COMPLETED`, `NO_EVIDENCE_FOUND`, `ACCESS_BLOCKED`, or `NOT_APPLICABLE` with a
reason. `COMPLETED` requires a real exact-schema table; the other states cannot
smuggle in rows. Output and provenance are immutable rather than overwritten.

`scripts/77_collate_interpretation.py` revalidates every task and creates:

- `results/atlas/regulatory_elements.tsv`;
- `results/atlas/cell_types.tsv`;
- `results/atlas/pathways.tsv`;
- `results/atlas/causal_tests.tsv`; and
- `results/tables/interpretation_coverage.tsv`.

The regulatory adapters uniquely lift GRCh37 variants for GRCh38 GENCODE,
SCREEN, and HOCOMOCO queries. Promoters may name only a supported same-locus
gene; SCREEN interval and HOCOMOCO motif rows remain unlinked
(`target_gene_id=NA`). The native-GRCh37 Nasser et al. ABC atlas is streamed
once into a checksum-bound overlap cache, excludes self-promoter links, uses
only 52 explicitly frozen non-transformed biosample labels across the four
domains, and retains a target only when its exact symbol maps to one supported
same-locus atlas gene. Its released score floor is 0.015 and the prespecified
primary tier is 0.02. Neither tier is experimental proof of regulation.
The native-GRCh37 Javierre et al. PCHi-C matrix is likewise streamed once.
Only CHiCAGO scores of at least 5 are retained, reciprocal bait-to-bait links
are deduplicated by their maximum score, and a promoter name must map exactly
to one supported same-locus atlas gene. Its 17 primary hematopoietic contexts
serve the immune domain; brain, metabolic, and vascular tasks are recorded as
`NOT_APPLICABLE`, not as negative evidence. A qualifying contact is evidence
of physical proximity, not proof of enhancer activity or causal regulation.

The FUMA single-cell strategy is frozen to the public
`vufuma/FUMA_scRNA_data` repository at commit
`dd526163ea80af1a80a6cdc80db167144500694b`. Exactly two matrices per domain
were selected before results were viewed. MAGMA v1.10 uses the FUMA Ensembl
v92 GRCh37 coding-gene boundaries (20,260 genes), a 1 kb upstream/downstream
window, the SNP-wise mean gene model, and one-sided positive gene-property
tests conditioned on each matrix's `Average` expression column. The bundled
MAGMA source was compiled natively for ARM64 with Eigen vectorization disabled;
the compiler identity and executable checksum are pinned. Vascular matrices
originate from mouse vascular studies already mapped by FUMA to human Ensembl
identifiers, so cross-species interpretation remains an explicit limitation.

The complementary chromatin-accessibility strategy is frozen to CATlas
Mendeley Data release v4 (`10.17632/yv4fzv6cnm.4`). Its three-file source
bundle contains 615,998 adult nuclei across 111 source cell types, 890,130
adult-present cCREs, and 111 cell-type-restricted peak sets. Exactly 43 adult
cell types were selected before enrichment results: 10 brain, 9 immune, 9
metabolic, and 15 vascular. For each trait, a one-sided hypergeometric test
compares genome-wide-significant variants with a fixed common, autosomal,
MAF-at-least-0.01, LD-pruned HapMap3 EUR universe, conditioning the background
on adult cCRE overlap. Each trait must cover at least 50% of that fixed universe,
matching the independently frozen FUMA reference-overlap floor; only variants
with valid trait P values enter its analyzable background, so absence is never
treated as a null association. The GRCh37 reference variants are uniquely
lifted to GRCh38 before interval overlap, zero-signal families are retained,
and BH is applied across all 43 cells per trait. This coarse annotation cannot
establish that a cell type, element, variant, or nearby gene is causal.

The current local production cache contains 72,394 uniquely lifted variants
after LD pruning, including 12,585 variants in the adult-cCRE background. All
43 cell types have nonzero fixed-reference overlap. Checksum-bound trait caches
are complete for 45/45 traits; coverage ranges from 64.7% to 99.8%. Two
coordinate-discordant telomere-length rsIDs were conservatively excluded and
recorded. Enrichment results remain intentionally uncomputed until the complete
interpretation task family can receive its pre-result lock.

Cell-type and pathway P values receive BH correction within their frozen test
families. Only FDR-supported rows enter those two canonical evidence tables.
The complete null family remains auditable through coverage. Promoter,
enhancer-promoter-link, and 3D-contact rows must link a fine-mapped same-locus
variant to an already supported gene; interval-only enhancer, open-chromatin,
and motif evidence may remain unlinked. Nearest-gene-only assignment is
forbidden. MR rows need at least three
instruments, minimum F ≥ 10, and all required diagnostics. Robust causal
wording additionally needs corrected multi-estimator agreement, clear
pleiotropy diagnostics, correct Steiger direction, characterized or adjusted
sample overlap, and compatible signal-level colocalization.

The causal runtime is frozen in
`config/interpretation_causal_runtime.json`. It uses R 4.3.3 aarch64 with an
isolated `.mr-env/library`, exact source archives for TwoSampleMR 0.7.9,
MR-PRESSO 1.0, cause 1.2.0, and lhcMR, 63 checksum-pinned dependency sources,
and an exact 153-package installed dependency closure. Stable PLINK
1.9.0-b.7.11 performs local-only LD operations against the same
checksum-pinned 503-sample 1000 Genomes Phase 3 EUR GRCh37 archive used by the
MAGMA workflow. The reference is materialized as a streamed SNP-major subset,
so the 3.60-GB archive never needs to be fully extracted. Palindromic variants
are always dropped (`harmonise_data` action 3), proxies and remote API clumping
are forbidden, and task seeds are deterministic. CAUSE is the primary
correlated-pleiotropy estimator when at least 100,000 pruned null-parameter
variants are available. lhcMR is sensitivity-only for material known overlap
and must receive real LDSC intercepts and real bidirectional MR starting
estimates; its random fallback paths are forbidden.

The public pathway adapters use the official Reactome v97 human GMT and the
archived GO 2026-08-05 release (go-basic ontology plus human UniProt GAF).
They also use the complete official human MSigDB v2026.1 symbols GMT and the
exact FUMA 1.5.6+ Ensembl-mapped MSigDB v2023.1Hs file distributed for MAGMA.
Source symbols map only when they resolve exactly and uniquely in GENCODE v26.
GO `NOT` annotations are excluded and remaining annotations are propagated only
through `is_a` and `part_of`. For each resource and locus, a one-sided
hypergeometric overrepresentation test compares high-confidence convergent
genes with the full eligible resource background. Every 10–1000-gene set is
retained in the unadjusted result, including zero-overlap null rows, so the BH
family is complete. This is pathway annotation, not evidence of pathway
activity, direction, mediation, or causality.

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
repository state, 18 of 20 interpretation source families are ready: GENCODE
promoters, SCREEN enhancer/open-chromatin layers, HOCOMOCO H14CORE, the ABC
2021 enhancer-gene atlas, the Javierre 2016 immune PCHi-C atlas, the pinned
FUMA single-cell/MAGMA bundle, the immutable CATlas adult scATAC bundle, the
Reactome v97 and GO 2026-08-05 pathway families, the official MSigDB v2026.1
human collection, the FUMA-prepared MSigDB v2023.1Hs MAGMA gene-set file, the
within-workflow regulatory cell-type layer, and all five causal estimator
families. LDSC-SEG remains blocked on its requester-pays official source, and
the molecular-QTL-derived cell strategy remains correctly blocked on upstream
molecular artifacts. The compressed
official 1000 Genomes Phase 3 GRCh37 European MAGMA
reference is checksum-pinned, but its 3.60 GB extracted members are deliberately
not materialized wholesale. Causal LD work uses a disk-bounded streamed subset;
the current local derivation retained 1,206,400 exact allele-matched HapMap3
variants in 232 MiB and passed a full PLINK allele-frequency read over all 503
reference samples. Five HapMap3 IDs had allele mismatches and 10,906 were absent
from the reference, all recorded in immutable derivation provenance.
Full MAGMA materialization remains guarded until at least 4,674,439,651 bytes
are free and requires an explicit large-extraction acknowledgement.
Four of six upstream canonical inputs do not yet exist. No downstream task
result has been generated or claimed.

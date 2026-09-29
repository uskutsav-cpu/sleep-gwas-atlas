# Brain6 prospective positional enrichment v2: full-set matching audit

**Decision: no full-set enrichment p-values.** The [prospectively frozen protocol](PROTOCOL.md) required at least 10 eligible, exact-width, chromosome-matched null anchors for *every one* of the 20 fixed geographic PLACO candidate regions. Two regions failed that pre-outcome gate (4 and 7 anchors). No pathway or tissue enrichment test was calculated, no caliper was relaxed, and no region was silently removed. This is an exploratory feasibility assessment; it does not change the frozen PLACO, LAVA, QTL, or manuscript classifications.

## Frozen set and source coverage

The input is exactly [25 pair-specific PLACO candidates](../loci/five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv) represented by [20 unioned geographic regions](../loci/five_track_cross_pair_region_reconciliation_v1/geographic_region_groups.tsv). The [source receipts](source_receipts.tsv) record SHA-256, byte size, and origin for these tables; all five complete pairwise PLACO tested-variant files; the 2,495 pinned GRCh37 LAVA blocks; [GENCODE v19 GRCh37 gene bodies](https://www.gencodegenes.org/human/release_19.html); [GTEx v8 median TPM across 54 bulk tissues](https://gtexportal.org/home/datasets/); and the [Reactome human pathway GMT](https://reactome.org/download-data). Public source files are held on the mounted Extreme SSD at `/Volumes/Extreme SSD/brain6-work/brain6-exploratory-enrichment-v2/inputs/`.

The complete PLACO files yielded **20,370,452 `TESTED` rows** across five tracks (5,514,399 insomnia–ADHD; 1,122,358 insomnia–MDD; 6,305,663 long sleep–bipolar; 1,131,670 long sleep–Parkinson; 6,296,362 long sleep–schizophrenia). **2,474/2,495** LAVA blocks contain at least one tested variant. GENCODE has **19,430** autosomal protein-coding gene bodies under the fixed parser; **19,153** overlap tested blocks. The 20 candidate regions overlap **177 distinct positional gene bodies**; **19/20** overlap at least one, and **19/20** have at least one GTEx-mapped gene. GTEx maps **18,288** of the 19,430 GENCODE gene IDs across 54 bulk tissues. Reactome has **2,868** source pathways, of which **1,680** have 10–500 uniquely symbol-mapped genes in the tested-block gene universe. These are input and mapping counts, not enrichment outcomes or causal assignments. [Region coverage](positional_annotation_coverage_20.tsv) and the [source-count summary](annotation_feasibility.json) preserve the details.

## Exact matching result

Each candidate region was matched to windows of its exact base-pair width, centered on the same chromosome's pinned LAVA blocks, excluding all 20 candidate overlaps. Prespecified matching calipers then required a 0.5–2 ratio for anchor LD-block span and every relevant PLACO-pair tested-variant count, plus the fixed protein-coding gene-count caliper. The [20-row matching audit](matching_feasibility.tsv) reports every candidate's covariates, surviving anchors, and sequential exclusion reason. The 18 regions that passed had **11–103** eligible anchors.

| Full-set failing region | Width | Positional genes | Same-chromosome anchors | Eligible | Sequential losses after edge/candidate exclusion |
|---|---:|---:|---:|---:|---|
| chr3:52,219,816–53,219,816 | 1,000,001 bp | 32 | 175 | **4** | LD-span 85, tested-density 1, gene-count 79 |
| chr17:43,049,526–45,343,136 | 2,293,611 bp | 29 | 65 | **7** | LD-span 1, tested-density 48, gene-count 1 |

The exclusion columns are sequential, mutually exclusive accounting within each target's ordered filter; they are not independent biological causes. The chr3 region's high positional gene count leaves few gene-count-matched controls after the LD-span filter. The chr17 region is unusually wide and its two relevant PLACO tracks leave few density-matched controls. Both remain members of the fixed 20-region universe. The full-20 test remains **NOT_ESTIMATED**.

## Modalities and interpretation

- **Pathways and bulk tissues:** Source and candidate coverage are adequate for annotation, but the required joint 20-region matched-null design fails. There are **zero** tested Reactome pathways and **zero** tested GTEx tissues in v2, hence no nominal or adjusted enrichment p-values.
- **Cell types:** The 54 GTEx columns are bulk tissues; **zero** source-versioned cell-type marker panels are available under this protocol. No cell-type enrichment is estimated.
- **Regulatory features:** The earlier exact-lead Ensembl source contains **27 lead queries and six coordinate-overlapping feature rows**; this branch has **zero matching regulatory queries for the 2,495 genome-wide null anchors**. It cannot provide a calibrated regulatory enrichment test. Exact-lead features remain descriptive only.
- **Genes:** Every gene label here means gene-body positional overlap under GENCODE v19. No colocalization, effector-gene, direction-of-effect, or causal interpretation follows.

The [machine-readable status](ENRICHMENT_STATUS.md), [preparation summary](preparation_summary.json), and [script](../../scripts/build_exploratory_enrichment_v2.py) allow the matching decision to be reproduced. The script reads only PLACO `CHR`, `BP`, and `status` for background matching; association P/Q columns are not consulted. The v2 protocol was committed as `cfb19bc1` before candidate gene mapping or matching results were calculated.

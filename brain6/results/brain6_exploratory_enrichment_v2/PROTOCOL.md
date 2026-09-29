# Prospective positional enrichment protocol (v2)

Frozen 2026-09-28 before computing candidate gene overlaps or enrichment outcomes. This is an exploratory follow-up of the fixed PLACO candidate set. It does not repair LAVA, identify causal genes, or qualify any confirmatory claim. The 25 pair-specific candidates are represented by the 20 unioned geographic intervals in `../loci/five_track_cross_pair_region_reconciliation_v1/geographic_region_groups.tsv`; a region is one observation even when several pairs contributed.

## Frozen sources and coordinates

- GRCh37/hg19 coordinates. Input region table above and the 2,495 LD-derived LAVA blocks in `ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile` are pinned repository inputs.
- Gene bodies: GENCODE release 19 `gencode.v19.annotation.gtf.gz`, reference chromosome annotation, `gene` records with `gene_type "protein_coding"` (or `gene_biotype` if applicable), on autosomes 1–22. Ensembl gene ID is the primary key; a region or null window maps to every gene body with at least one base of overlap. No distance-to-lead, eQTL, or causal prioritization is inferred.
- Tissue panel: GTEx Analysis v8 `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`, all 54 bulk tissue columns. Match Ensembl gene IDs after removing version suffix; missing genes are omitted and counted. Each gene's tissue score is `log2(1 + median TPM in tissue) - mean_across_54_tissues(log2(1 + median TPM))`.
- Pathways: official Reactome human `ReactomePathways.gmt.zip`, release file dated 2026-06-21; map gene symbols only when they map uniquely to a GENCODE v19 protein-coding gene ID. Use all pathways with 10–500 mapped genes in the genome-wide tested-window gene universe. Nested pathways remain distinct tests.
- PLACO tested background: the complete five `brain6/results/placo/<pair>/variants.tsv.gz` files, counting only rows with `status=TESTED` and a valid chromosome and BP. Candidate-specific matching uses each pair listed for that geographic region; no P or Q values are read for matching.

Source URLs, SHA-256, file sizes, and exact parser counts will be recorded in `source_receipts.tsv`. Downloaded public files remain on the SSD at `/Volumes/Extreme SSD/brain6-work/brain6-exploratory-enrichment-v2/inputs/`.

## Background and eligibility

For each of the 20 target intervals, construct one null interval of *exactly the same base-pair width* centered at each LAVA block midpoint on the *same chromosome*. A null is eligible if it lies within the chromosome extent spanned by the LAVA blocks and does not overlap any target region. The target's LD-block proxy is the length of the LAVA block containing its midpoint; each null anchor's block length must be within 0.5–2 times that length. For every contributing PLACO pair, the null's `TESTED` variant count must be within 0.5–2 times the target count; zero requires zero. Its positional protein-coding gene count must differ by no more than `max(2, ceil(0.5 * target_count))`; for a zero-gene target require zero genes. These are prespecified calipers, not adaptive to enrichment results. Require at least 10 eligible distinct null anchors for every target; otherwise no pathway or tissue p-values are issued and the shortfall is reported per region.

Generate 10,000 accepted genome-wide null replicates using seed 20260928. Each replicate samples one eligible null independently for every target, rejecting and resampling if two null intervals overlap on the same chromosome. Use up to 200,000 attempts for 10,000 accepted replicates; if that cannot be achieved, issue no enrichment p-values and report the acceptance count. A target whose midpoint is outside all pinned LAVA blocks also stops inference. The eligible pools, matching covariates, and all exclusion reasons are output before statistical calculations and retained.

## Tests and reporting

- **Bulk tissue expression:** For each interval, average each of the 54 gene tissue scores across its overlapping GTEx-mapped genes (zero for a gene-empty interval). The statistic is the sum of these 20 interval means. Compare each tissue against the 10,000 matched-null statistics, one-sided high-tail empirical `p=(1+#null>=observed)/(1+B)`. Benjamini–Hochberg correction across all 54 tissues. Report all tissues, nominal and adjusted p, mean score, and mapping coverage.
- **Reactome pathways:** For each pathway, count the number of the 20 intervals with at least one overlapping mapped pathway gene. Compare to the analogous null count, with the same empirical p formula and Benjamini–Hochberg across *all* pathways meeting the 10–500 universe-gene criterion. Report every tested pathway, even if observed count is zero, plus mapped-gene and region coverage.
- **Regulatory and cell-type:** No enrichment p-value from these panels. The available lead-variant Ensembl regulatory queries lack identical genome-wide null queries, and GTEx median TPM represents bulk tissues, not cells. Quantify the available annotations and untested background in the audit. Neither descriptive annotations nor bulk tissue scores are called cell-type enrichment.

No threshold, mapping rule, source, or set will be changed after looking at output. If matching fails, report a numerical feasibility audit with the fixed protocol intact.

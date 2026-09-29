# Brain6 exploratory functional follow-up: Phases 11–14

**Status: source-bound exploratory results only; no candidate promotion.** The universe is exactly 25 pair-specific PLACO candidates in 20 geographic regions. Canonical LAVA remains `FAILED_QC_NOT_PROMOTED`. This branch does not supply a final shared-locus tier, a causal gene, or independent two-trait locus replication. The [protocol](PROTOCOL.md) fixed the source, N, region, QTL, and coloc gates before new regional QTL association values or posterior results were inspected.

## Source and method accounting

- Four [official eQTL Catalogue release-7](https://www.ebi.ac.uk/eqtl/Data_access/) GTEx brain panels were inspected: QTD000171 cortex eQTL, QTD000175 cortex Leafcutter sQTL, QTD000176 frontal-cortex eQTL, and QTD000180 frontal-cortex Leafcutter sQTL. The [Catalogue license](https://www.ebi.ac.uk/eqtl/License/) is CC BY 4.0. Its [column documentation](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md) defines GRCh38 REF/ALT, ALT-effect beta, SE, P, MAF, and autosomal analyzed N as `AN/2`. Twenty-four tabix regional slices and receipts live on the SSD under `/Volumes/Extreme SSD/brain6-work/brain6-exploratory-functional-v1/qtl_regional/`; each receipt records official URL, ETag, Last-Modified, source size, GRCh37 and GRCh38 intervals, row count, and local SHA-256. These are **regional** hashes, not hashes of the four full source files. All four previous published molecular credible-set files and lead lookups were reused without rewriting.
- The exact Jansen insomnia and Demontis ADHD **raw** releases were checked against their prior source hashes. They supplied literal per-SNP `N` and per-SNP `Nca`/`Nco`, respectively, for twelve new regional slices. The original Brain6 harmonized `N` fields are effective-N transforms and were not used. Source-level case fraction is approximate for the insomnia variants, while ADHD variant-specific case and control counts are present. The resulting GWAS–QTL calculation is a component-GWAS sensitivity, independent of canonical LAVA admission.
- Other 19 pair-specific candidates have unresolved original-GWAS N/model or cohort-mixture requirements for this new coloc analysis. These were method-held before inspecting their QTL association values. The six insomnia–ADHD candidates are fully represented in the [264-row gate table](coloc_gate.tsv): four eligible tests, 169 without a strong shared-variant QTL signal, 71 with fewer than 500 uniquely shared variants, four without a strong shared-variant GWAS signal, and 16 panel/window combinations without QTL records.

## GWAS–QTL coloc sensitivity

Only the **ADHD component** of the `chr5:87185500–88185500` insomnia–ADHD candidate met every frozen gate. Four single-signal `coloc.abf` 5.2.3 tests used 1,505–1,624 exact GRCh38 position plus unordered-allele matched variants. Priors were `p1=10⁻⁴`, `p2=10⁻⁴`, `p12=10⁻⁵`; QTL `sdY` was estimated by coloc from MAF, SE, and `AN/2` because a direct molecular phenotype SD was unavailable. The [result table](coloc_results.tsv) retains all H0–H4 posteriors and source input paths:

| Molecular context | Molecular trait | Shared variants | PP.H3 | PP.H4 |
|---|---|---:|---:|---:|
| Cortex eQTL | ENSG00000271904 | 1,624 | 0.999998 | 0.00000218 |
| Cortex sQTL | `5:88270585:88436273:clu_38358_+` | 1,624 | 0.99999996 | 0.0000000316 |
| Frontal-cortex eQTL | ENSG00000250377 | 1,505 | 0.999863 | 0.00000939 |
| Frontal-cortex eQTL | ENSG00000271904 | 1,505 | 0.785945 | 0.082385 |

The strongest PP.H4 is **0.0824**, so none of these four single-signal tests provides persuasive support for an ADHD–QTL shared variant under the specified priors. This is **not** a negative result for the insomnia–ADHD pair: no insomnia-component coloc passed the predeclared molecular and GWAS signal gates, and the two-trait locus was not independently replicated. One-causal-variant ABF is especially fragile in complex loci with multiple association signals; the QTL `sdY` approximation and ADHD meta-analysis sample heterogeneity add uncertainty. These posterior values are limited exploratory sensitivities, not calibrated evidence against all regulatory mechanisms. The molecular trait IDs are not assigned causal-gene status.

## Replication, regulatory and enrichment status

The [25-candidate table](candidate_functional_status_25.tsv) and [20-region table](geographic_region_status_20.tsv) carry every original candidate and geographic group. The prior published QTL credible-set memberships and exact-lead Ensembl regulatory overlaps are descriptive only; pair-counted memberships recur across shared geographic regions and cannot serve as independent enrichment hits. Six exact-lead regulatory feature rows were retained from the prior official Ensembl snapshot. An empty exact-lead query does not exclude regulatory activity elsewhere in its region.

The archived FinnGen R13 ADHD genetic-correlation result is **global pair-level directional** context; it cannot establish locus-specific two-trait replication. Existing candidate long-sleep adjacent-phenotype lookups use ≥10 h, UKB-overlapping data and also do not meet this gate. No additional verified, compatible, cohort-disjoint regional replication input was acquired. Thus all 25 rows retain `independent_two_trait_locus_replication=NOT_ESTABLISHED`.

Tissue/cell-type and pathway enrichment remain `NOT_ESTIMATED`: no source-reviewed causal/weighted gene set, tested gene universe, LD/region-aware background, tissue or cell panel, ontology collection, and multiplicity family are jointly specified. Published QTL credible-set membership, exact-lead coordinate features, and nearest genes cannot be relabeled as such enrichment. The [source gate](source_gate.tsv) and [provenance](provenance.json) retain explicit missing-input states and hashes.

No canonical output, source archive, prior exploratory snapshot, replication classification, or protected candidate slot was modified. No outcome was used to choose a QTL panel, priors, gene, candidate region, or method gate.

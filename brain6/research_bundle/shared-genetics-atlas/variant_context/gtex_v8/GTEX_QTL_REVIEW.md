# GTEx v8 measured QTL context for Brain6 representatives

Queried 27 October 2026 using the [GTEx Portal API](https://gtexportal.org/api/v2/redoc), exact GRCh38 REF/ALT alleles and the 20 variant–region representatives listed in `../selected_variant_gene_context.tsv`. For each variant, the significant single-tissue eQTL and sQTL endpoints were queried across all returned tissues. All 40 requests succeeded; the API returned 76 rows (75 unique complete association records; the response for rs324886 includes one duplicated adipose sQTL record). Fourteen of the 20 variants have at least one returned row. The stored table retains the duplicate exactly rather than silently changing the source. Empty results mean no record was returned by that endpoint for this variant and release, not that the association was tested with adequate power or has no effect. This is a selected, post-hoc lookup, not a systematic test for enrichment or independent confirmation of a shared GWAS causal variant.

For eQTLs, GTEx identifies significant variant–gene pairs using gene-specific thresholds tied to its eGene discovery procedure; the displayed `pValue` is nominal, not a per-row BH-adjusted q-value. See [GTEx methods](https://gtexportal.org/home/methods). Its normalized effect size (NES) is an ALT-versus-REF regression slope in normalized expression space, not a biological fold change. sQTL NES relates to the stated splice phenotype; it should not be interpreted as gene expression or direct activity. Tissue, gene and source version matter.

## Exact rs62264764 results

The A alternate allele at **chr3:117920728 G>A** has a GTEx v8 significant-pair eQTL record in **adrenal gland** for *LSAMP* (ENSG00000185565.11; nominal p=6.86×10⁻⁷, NES=+0.483) and RP11-384F7.2 (ENSG00000239268.2; nominal p=1.09×10⁻⁵, NES=+0.424). The latter shares its unversioned Ensembl gene ID with the current annotation LINC03051. These are measured variant–expression associations in bulk adrenal tissue, not proof that either gene mediates sleep or ADHD. Differences between adrenal and brain tissue are central to interpretation. The single-tissue significant sQTL endpoint returned no rows for this variant.

To avoid conflating rs62264764 with nearby rs62264767, we separately used GTEx's dynamic eQTL endpoint for the exact allele and both annotated genes in cortex, frontal cortex and adrenal gland. This service returns nominal tests regardless of significance; the table below contains its output. The small group sizes, post-hoc selection and multiple gene/tissue comparisons limit interpretation.

| Gene (GTEx v8) | Tissue | n | ALT NES | Nominal p |
|---|---|---:|---:|---:|
| LSAMP | Cortex | 205 | −0.0076 | 0.870 |
| LSAMP | Frontal cortex BA9 | 175 | +0.0503 | 0.355 |
| RP11-384F7.2 | Cortex | 205 | +0.1361 | 0.0954 |
| RP11-384F7.2 | Frontal cortex BA9 | 175 | −0.1658 | 0.0323 |
| LSAMP | Adrenal gland | 233 | +0.4828 | 6.86×10⁻⁷ |
| RP11-384F7.2 | Adrenal gland | 233 | +0.4238 | 1.09×10⁻⁵ |

The frontal RP11-384F7.2 p=0.0323 is unadjusted and was not returned in the significant-pair query. It is a weak exploratory lead, not an established brain eQTL; multiplying by six post-hoc queries gives 0.194 as a descriptive Bonferroni value, but this was not a prospectively frozen six-test family. The other brain tests do not establish absence of effects. Dynamic adrenal results match the corresponding significant-pair records within rounding. Full summaries, counts, exact URLs and response hashes are in `rs62264764_dynamic_summary.tsv`. The dynamic API also returned per-sample arrays; these were not saved, and the hash records the response as received, so exact response bytes are not available in this package.

## Other selected candidates

- **rs4916723**, chr5:88558577 A>C: the significant-pair endpoint returns a cortex eQTL for TMEM161B-AS1 (nominal p=2.50×10⁻⁸, ALT NES=−0.407), along with other brain regions. It is an ADHD-only member of the selected credible-set union in this pilot, with very low insomnia PIP. Thus this is molecular context for an ADHD-prioritized marker, not evidence of a shared insomnia–ADHD causal variant. rs324886 has a cortex TMEM161B-AS1 record too (p=1.25×10⁻⁵, NES −0.330). Both also have selected sQTL records in non-brain tissues.
- **rs1476535**, chr7:114430980 C>T, which overlaps FOXP2, has significant-pair eQTL records for FOXP2 in pituitary, adrenal and subcutaneous adipose tissues; no brain record returned. Its published shared-locus association remains distinct from a functional or causal gene assignment.
- **rs2431108 and rs77960**, chr5 at the other prioritized region, have eQTL and sQTL records for RP11-6N13.1 in testis, but no significant brain record returned in these queries. Existing lead-level brain nominal results and the present tissue records do not constitute GWAS–QTL coloc.

The `query_status.tsv` gives all 40 exact query URLs, response counts, and SHA-256 hashes; each static source response is preserved under `responses/`. `significant_associations.tsv` reproduces all returned rows. The source-level GTEx thresholds do not correct the multiple model- and phenotype-driven follow-up decisions made in this project. None of these variant-level QTL records establishes which correlated regional variant drives a molecular trait, nor whether it shares the GWAS causal signal. Proper follow-up requires dense per-gene QTL summary statistics, ancestry-matched LD, and multi-signal GWAS–QTL colocalization in relevant tissue, with independent experimental validation.

The later [dense rs62264764-region eQTL–GWAS comparison](../rs62264764_dense_qtl/RESULTS.md) finds low shared-signal posterior support for the two selected genes and tissues under a single-signal model, despite the adrenal exact-variant associations. This is sensitive to coverage and assumptions and does not rule out secondary sharing.

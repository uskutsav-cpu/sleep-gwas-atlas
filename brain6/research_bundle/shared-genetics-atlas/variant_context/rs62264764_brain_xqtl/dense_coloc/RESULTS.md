# rs62264764 region: oligodendrocyte expression and GWAS

In a **post-selected single-signal colocalization model**, the strongest oligodendrocyte LINC03051 expression association in the analyzed chromosome 3 window appears distinct from the dominant insomnia and ADHD associations. The model favors H3 (both traits associated, different dominant signals) over H4 (shared dominant signal) for each GWAS.

| Comparison | Matched SNPs | H3, p12=10⁻⁵ | H4, p12=10⁻⁵ | H4, p12=10⁻⁶ |
|---|---:|---:|---:|---:|
| LINC03051 expression–insomnia | 2,458 | 0.9955 | 0.00167 | 0.000167 |
| LINC03051 expression–ADHD | 2,340 | 0.9813 | 0.00438 | 0.000440 |

These are model-conditional posteriors, not corrected p-values or probabilities that a gene mediates disease. The analysis does **not** rule out a secondary shared signal: both gene expression and disease GWAS can contain multiple signals, and the model assumes at most one causal variant per trait. It also does not contradict the measured rs62264764 eQTL: a marker can tag an expression signal without being the shared causal variant. An exact-variant association alone does not establish mediation or validate the AlphaGenome DNase prediction. No AlphaGenome values entered this model.

The strongest matched expression association in both fits is rs4687889 (chr3:117818592 T>C, GRCh38; nominal p=2.28×10⁻²²); the best matched GWAS markers differ. The first two markers in the conditional H4 ranking are rs4687889 and rs4688051 for the insomnia fit, but conditional SNP weights describe the model **given H4**, not absolute causal probabilities. These are descriptive observations from selected data, not independent replication or proof of distinct causal variants.

## Data and validation

The combined ROSMAP oligodendrocyte data are from aging-brain single-nucleus expression, [NIAGADS NG00184.v1](https://dss.niagads.org/datasets/NG00184/), not an ADHD or sleep cohort. The chromosome 3 SNP file and index were retrieved from the uncompressed publisher archive using the documented [byte ranges](archive_ranges.tsv) and matched publisher MD5 checksums; the full archive was not downloaded or hashed. The index has an older timestamp than the data file because of retrieval order, but its bytes matched the publisher's checksum and indexed reading worked. The exact rs62264764 beta, SE and p matched the earlier portal record.

Within the selected window, 4,020 gene–variant records were found; 3,402 had valid nonpalindromic statistics. Of those, 2,458 matched insomnia and 2,340 matched ADHD on unique lifted position and unordered alleles. The original GWAS files lack allele frequencies, so frequency compatibility was not independently checked. The one-megabase window may truncate cis/TAD effects, and excluded or missing variants can affect posteriors. The [matching receipt](matching_receipt.json), [protocol](PROTOCOL.md), [preparation code](prepare.py), [full posterior table](posteriors.tsv), [conditional ranking](top_conditional.tsv), [R model](run.R), and [independent arithmetic verification](verification.json) record the run. All four posterior vectors were reproduced to within 4.5×10⁻¹⁶; this checks arithmetic, not calibration or assumptions. Source files and restricted GWAS-derived matched records remain on the SSD.

This is an oligodendrocyte test only. The earlier portal also reported an excitatory-neuron eQTL; the full excitatory file was not acquired because the other SSD image could not create files. A later [targeted regional retrieval and excitatory pilot](../excitatory_coloc/RESULTS.md) avoids that full-file download. Study overlap, age, ancestry, effect-prior choice, incomplete support and LD–Z consistency remain important. Appropriate follow-up requires multi-signal methods with study-matched LD and QC, independent disease cohorts and experiments. The failed LAVA gates remain failed.

Data for this study were prepared, archived, and distributed by the National Institute on Aging Alzheimer's Disease Data Storage Site (NIAGADS) at the University of Pennsylvania (U24-AG041689), funded by the National Institute on Aging. We acknowledge the ADSP FunGen xQTL Atlas contributors and research participants; accession **NG00184**.

A later [brain-cell LD–Z readiness diagnostic](../ldz_readiness/RESULTS.md) found notable residuals despite similar allele frequencies in the reference; a validated multi-signal fit remains outstanding.

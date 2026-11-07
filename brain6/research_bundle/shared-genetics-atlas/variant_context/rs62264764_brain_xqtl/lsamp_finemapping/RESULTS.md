# Study-produced LSAMP fine mapping and brain-cell QTL follow-up

The [NIAGADS ADSP FunGen xQTL Atlas](https://dss.niagads.org/datasets/NG00184/) reports **LSAMP expression credible-set records** in both oligodendrocytes and excitatory neurons in the combined ROSMAP aging-brain dataset. This adds brain-cell-specific molecular evidence for LSAMP in the rs62264764 region. It does not establish LSAMP as the gene mediating insomnia or ADHD.

| Cell type | Reported LSAMP 95% credible-set records | Highest-PIP marker and PIP | rs62264764 in reported set? |
|---|---:|---|---|
| Oligodendrocytes | 10 | chr3:117871355 C>T, 0.271 | No |
| Excitatory neurons | 3 | chr3:117793105 A>G, 0.499 | No |

Coordinates are GRCh38. These are the publisher's **single-context** fine-mapping records, filtered to HMT-significant variant–target pairs with additional set QC; they are not necessarily every record from the original fit. PIPs and credible sets depend on the publisher's model and LD assumptions and are not experimental proof. Their particular LD construction has not been independently verified here. Neither rs62264764 nor its strongly linked neighbor rs62264767 appears in these two reported sets; that does not exclude another signal or establish either variant's causal status. No LINC03051 *target* records were found in the four inspected chromosome 3/cell-type CS95 files, which is not evidence that it has no expression association or that all fine-mapping elsewhere is negative. Some LSAMP rows annotate LINC03051 as the **variant-mapped gene**; that annotation is not the molecular target.

## Expression–GWAS comparison

Using all tested regional QTL records rather than only significant records, a separate **single-signal** ABF model gave:

| Cell and GWAS | Matched SNPs | H3, distinct, p12=10⁻⁵ | H4, shared, p12=10⁻⁵ | H4, p12=10⁻⁶ |
|---|---:|---:|---:|---:|
| Oligodendrocyte LSAMP–insomnia | 2,458 | 0.9961 | 0.00106 | 0.000106 |
| Oligodendrocyte LSAMP–ADHD | 2,340 | 0.9731 | 0.0127 | 0.00129 |
| Excitatory LSAMP–insomnia | 2,458 | 0.9971 | 0.00000865 | 0.000000865 |
| Excitatory LSAMP–ADHD | 2,340 | 0.9856 | 0.0000681 | 0.00000681 |

H3/H4 are posteriors conditional on a selected region, priors, matched variants and a one-signal model; they are not corrected p-values or gene-mediation probabilities. They favor distinct dominant expression and disease signals under these settings but **cannot exclude shared secondary signals**. The two cells share donors, and four pairwise analyses do not establish a three-trait relationship. The exact rs62264764 LSAMP expression associations in the unfiltered files are nominal p=0.00224 in oligodendrocytes and p=0.0000601 in excitatory neurons; neither was returned as an HMT-significant exact-variant record in the prior portal lookup. The opposite signs are not a formal cell-type interaction test.

The study-produced oligodendrocyte LSAMP set shares three exact variant IDs with the prior GTEx adrenal LSAMP SuSiE credible set; the excitatory set shares none. These are post-hoc set intersections, **not** colocalization. That earlier [adrenal multi-signal pilot](../../rs62264764_dense_qtl/multisignal_pilot/RESULTS.md) was not validated because of LD–Z and other source limitations. The present ABF models and study-produced QTL sets do not rescue it. Prior [brain-cell LD diagnostics](../ldz_readiness/RESULTS.md) also found notable residuals with the small out-of-sample reference.

## Provenance and limits

The publisher's complete CS95 archive and manifest matched their MD5 values; the four selected chromosome 3 cell/SNP-or-indel members also matched their published MD5s. The [extracted target rows](study_credible_set_rows.tsv), [archive receipt](fine_mapping_receipt.json), [set verification](set_verification.json), [protocol](PROTOCOL.md), [matching receipts](matching_receipt_oli.json), [excitatory receipt](matching_receipt_exc.json), [posterior table](posteriors.tsv), [model code](run.R), and [arithmetic verification](verification.json) preserve the analysis. Eight posterior vectors were independently reproduced within 4.5×10⁻¹⁶; this is an arithmetic check, not statistical or biological validation. Full oligodendrocyte QTL files and selected excitatory regional bytes were previously checked as described in the [oligodendrocyte](../dense_coloc/RESULTS.md) and [excitatory](../excitatory_coloc/RESULTS.md) reports; the full excitatory chromosome source was not downloaded or whole-file hashed. The GWAS inputs lack allele frequencies, the one-megabase window can truncate cis/TAD signals, the QTL and disease sample overlap is unknown, and the data are from aging brain. The source files and restricted matched GWAS data remain on the SSD.

These findings do not validate an AlphaGenome prediction, establish a sleep–ADHD effector gene, demonstrate independent replication, or repair failed LAVA QC. Multi-signal tests require suitable validated LD, input-model checks and independent evidence, followed by experimental perturbation.

Data for this study were prepared, archived, and distributed by the National Institute on Aging Alzheimer's Disease Data Storage Site (NIAGADS) at the University of Pennsylvania (U24-AG041689), funded by the National Institute on Aging. We acknowledge ADSP FunGen xQTL Atlas contributors and participants; accession **NG00184**.

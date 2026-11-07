# rs62264764 region: excitatory-neuron expression and GWAS

A post-selected **single-signal** model favors distinct dominant LINC03051 expression and disease signals in the rs62264764 region, as it did in the earlier oligodendrocyte pilot. The results concern excitatory-neuron expression from aging-brain donors, not an ADHD or sleep cohort.

| Comparison | Matched SNPs | H3, distinct, p12=10⁻⁵ | H4, shared, p12=10⁻⁵ | H4, p12=10⁻⁶ |
|---|---:|---:|---:|---:|
| LINC03051 expression–insomnia | 2,458 | 0.9971 | 0.00000884 | 0.000000884 |
| LINC03051 expression–ADHD | 2,340 | 0.9856 | 0.0000699 | 0.00000699 |

H3/H4 are posterior probabilities conditional on the selected variants, effect priors and one-signal model, not adjusted p-values or evidence that a specific SNP or gene causes disease. Both traits may have multiple signals, so low H4 does not exclude secondary sharing. A variant-level eQTL at rs62264764 can tag a molecular signal even if the model's strongest disease and expression signals differ. This is not proof of distinct causal variants, independent replication, validation of AlphaGenome's DNase prediction, or resolution of the failed LAVA quality gates. AlphaGenome scores were not used in the model.

For context, the smallest matched nominal expression p-value is 1.54×10⁻²⁸ at rs9829042 (GRCh38 chr3:117956015 A>G) in both fits. It is a marker-level association, not a causal assignment. The [support comparison](support_comparison.json) confirms that the excitatory and prior oligodendrocyte fits use identical variant keys and GWAS effects within each trait. Their expression measurements share ROSMAP donors and are **not independent confirmations**. Opposite exact-variant expression directions between cell types are not a formal cell-type interaction test.

## Regional retrieval and validation

The [NIAGADS ADSP FunGen xQTL Atlas NG00184.v1](https://dss.niagads.org/datasets/NG00184/) publishes all tested associations as large tar archives. We retrieved the publisher's excitatory-neuron chromosome 3 tabix index, verified its MD5, and used its intervals and HTTP byte ranges to extract **1,762,845 compressed bytes** from the 671,965,359-byte chromosome file. The archive member's header and index matched the intended file, and the exact rs62264764 effect, SE and p matched the previously downloaded portal record. The extraction algorithm was separately checked against the complete checksum-verified oligodendrocyte file: every regional record matched pysam/tabix. Full excitatory chromosome bytes were not downloaded or hashed, and the selected range has no publisher-provided checksum; its SHA-256 records what was retrieved, not a whole-source authenticity check.

The region contained 4,020 target records, 3,402 of which passed valid nonpalindromic SNP/statistic checks; 2,458 and 2,340 matched insomnia and ADHD, respectively. The regional window may truncate cis/TAD signals, and GWAS source allele frequencies are unavailable. Population matching, GWAS–QTL overlap, multiple causal variants, effect priors and aging-brain relevance remain unresolved. A multi-signal analysis requires suitable LD and LD–Z diagnostics, and gene mediation needs independent evidence and experiments.

[Protocol](PROTOCOL.md), [range extraction and control](tbi_extract.py), [retrieval receipt](retrieval_receipt.json), [matching code](prepare.py), [matching receipt](matching_receipt.json), [posterior table](posteriors.tsv), [R model](run.R), [independent arithmetic verification](verification.json), and [support-check code](compare_support.py) document the work. Python reproduced all four posterior vectors within 4.5×10⁻¹⁶; this tests arithmetic, not model calibration. Source and restricted matched records are on the SSD. See the [oligodendrocyte pilot](../dense_coloc/RESULTS.md) for its distinct, fully checksum-verified source file.

Data for this study were prepared, archived, and distributed by the National Institute on Aging Alzheimer's Disease Data Storage Site (NIAGADS) at the University of Pennsylvania (U24-AG041689), funded by the National Institute on Aging. We acknowledge ADSP FunGen xQTL Atlas contributors and participants; accession **NG00184**.

A later [brain-cell LD–Z readiness diagnostic](../ldz_readiness/RESULTS.md) found notable residuals despite similar allele frequencies in the reference; a validated multi-signal fit remains outstanding.

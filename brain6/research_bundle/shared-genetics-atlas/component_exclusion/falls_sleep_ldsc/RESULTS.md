# Falls and insomnia: public GWAS sensitivity

A public [Neale Lab GWAS of UK Biobank field 2296](https://nealelab.github.io/UKBB_ldsc/h2_summary_2296.html) provides an additional physical-outcome sensitivity: it asks about falls in the last year and does **not** include a sleep or exhaustion question. Falls are influenced by many physical, neurological, environmental and other factors; they are not a frailty score and are not FFS with exhaustion removed.

| Insomnia source | Falls–insomnia genetic correlation r_g (SE) | LDSC p-value | Two-test Bonferroni p |
| --- | ---: | ---: | ---: |
| Brain6 questionnaire insomnia | 0.4302 (0.0346) | 1.50×10⁻³⁵ | 3.01×10⁻³⁵ |
| FinnGen registry insomnia | 0.2458 (0.0414) | 2.93×10⁻⁹ | 5.85×10⁻⁹ |

Both correlations meet this limited two-test benchmark; it does not correct for the wider atlas or earlier source searches. UK Biobank's [field coding](https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=100539) orders answers from no falls to one and then multiple falls, and the Neale page describes an ordinal phenotype. The exact source transformation and sign coding have not been independently reconstructed, so interpret the positive direction subject to that limitation. The published [Trajanoska et al. falls GWAS](https://doi.org/10.1038/s42003-020-01256-x) had already reported a positive falls–insomnia genetic correlation; this is not a novel discovery or an independent replication of the UKB-based signal. The Neale ordinal GWAS and that paper's binary falls GWAS are different analyses and should not be treated as identical.

The two calculations used different SNP supports and no paired statistical contrast was run; the difference between estimates cannot be attributed to a particular phenotype or cohort difference. The questionnaire insomnia and falls sources include UK Biobank participants, whereas the registry source is FinnGen; exact participant/relative overlap and the effects of ascertainment and ancestry are not verified. The unconstrained cross-trait LDSC intercepts were 0.0756 (SE 0.0062) and −0.0105 (0.0056), respectively; an intercept is not a participant-overlap count. Registry insomnia's univariate intercept was 1.0778 (0.0096), and the general-European LD reference and binary-trait N treatment may not fit FinnGen perfectly. The two correlations do not establish causality, local sharing, a sleep-independent physical-frailty construct or any sleep/exhaustion-excluded composite GWAS.

## Inputs and reproducibility

The falls file has 1,096,648 unique SNPs and constant N=360,344. The questionnaire fit retained 1,020,637 valid-allele SNPs, and the registry fit retained 1,053,614. A separate allele audit found no mismatched or palindromic allele pairs on the shared sets, and all three files passed complete gzip reads with no duplicate IDs or invalid rows under the scripted checks. The falls file's SHA-256 is recorded for identification; no independent publisher checksum was available. The questionnaire N=313,750.036 and registry N=497,916 are values in their processed files, not verified numbers of independent overlapping individuals or per-variant effective N.

[Protocol](PROTOCOL.md), [download and receipt](download_receipt.json), [input audit](input_audit.json), [audit code](audit.py), [LDSC runner](run.sh), and full [questionnaire](questionnaire.log) and [registry](registry.log) logs preserve the details. The compressed falls file remains on the SSD at `/Volumes/Brain6NHANES/neale-falls-20261028/falls.ldsc.bgz`. These are genome-wide correlations; no original LAVA gate was repaired and no locus or gene is established by them.

A subsequent [paired common-marker comparison](paired_source_sensitivity/RESULTS.md) estimates the difference between the questionnaire and registry correlations, subject to cohort, phenotype and LD-reference limitations.

# rs62264764: dense GTEx expression–GWAS comparison

The exact rs62264764 A allele has measured adrenal eQTL associations with *LSAMP* and RP11-384F7.2/LINC03051. However, in this post-selected **single-signal** regional model, neither adrenal expression signal shows strong shared-variant support with the Brain6 insomnia or ADHD GWAS. Frontal-cortex comparisons also show low shared-signal support, but are less informative because the model favors a GWAS signal without a detectable eQTL signal. This does not disprove secondary sharing or show that any particular gene is irrelevant.

| GTEx tissue | Expression trait | Brain6 GWAS | Matched variants | H2: GWAS only | H3: both, distinct | H4: shared | H4 with stricter prior |
|---|---|---|---:|---:|---:|---:|---:|
| Adrenal | LSAMP | insomnia | 1,632 | ~0 | 0.9961 | 0.00095 | 0.000095 |
| Adrenal | LSAMP | ADHD | 1,556 | ~0 | 0.9640 | 0.0174 | 0.00177 |
| Adrenal | RP11-384F7.2 | insomnia | 2,386 | ~0 | 0.9961 | 0.00095 | 0.000095 |
| Adrenal | RP11-384F7.2 | ADHD | 2,287 | ~0 | 0.9726 | 0.0132 | 0.00134 |
| Frontal cortex BA9 | LSAMP | insomnia | 1,632 | 0.8598 | 0.0968 | 0.0406 | 0.00422 |
| Frontal cortex BA9 | LSAMP | ADHD | 1,556 | 0.8510 | 0.0901 | 0.0407 | 0.00422 |
| Frontal cortex BA9 | RP11-384F7.2 | insomnia | 2,386 | 0.6946 | 0.1408 | 0.1622 | 0.0190 |
| Frontal cortex BA9 | RP11-384F7.2 | ADHD | 2,287 | 0.7240 | 0.1404 | 0.1230 | 0.0138 |

The H4 column uses p1=p2=10⁻⁴ and p12=10⁻⁵; the final column uses p12=10⁻⁶. The default eQTL effect prior SD is 0.15 and the case–control GWAS prior SD is 0.20. H2, H3 and H4 are model-conditional probabilities for GWAS-only, separate signals, and a shared signal, respectively, with the eQTL assigned trait 1. They are not empirical probabilities of a gene mediating disease. Across pre-specified alternative effect priors (SD 0.10 and 0.30 for both traits), none of the eight H4 values exceeded 0.163 over the tested shared priors. The 48 complete fits are in [posteriors.tsv](posteriors.tsv).

For context, rs62264764 itself has adrenal *LSAMP* nominal p=6.86×10⁻⁷ and ALT normalized expression slope +0.483; the strongest matched *LSAMP* adrenal eQTL in this interval has nominal p=3.94×10⁻³⁴. Thus the variant-level eQTL is not necessarily the strongest molecular signal. For *LSAMP* in frontal cortex, the exact variant's nominal p=0.355. For RP11-384F7.2 frontal cortex, its nominal p=0.0323 is post hoc, was not returned as a significant GTEx pair, and is not independent gene-level or tissue-specific evidence. The four exact-variant eQTL betas, SEs and p-values match the previously saved GTEx dynamic results within rounding.

## Interpretation and boundaries

The GRCh38 1-Mb window was derived from the previously chosen Brain6 GRCh37 region and can truncate gene cis windows. Depending on gene and trait, the matched variants cover a subset of eligible QTL rows; see [matching receipt](match_receipt.json). The frozen GWAS region files contain no allele frequencies, preventing an independent frequency-based harmonization check. Palindromic, invalid and ambiguous variants were excluded rather than treated as negative. GWAS p-values in the matched tables were reconstructed from beta/SE and may differ from source p-values. QTL rows with zero-frequency or invalid SE were retained in extraction and excluded from fits. Single-signal coloc cannot resolve multiple independent eQTL or GWAS signals. The overlapping genes and tissues are not independent replications, and the region was chosen after previous findings. Sample overlap, population differences, GTEx bulk-tissue limitations and uncertain effect priors also matter. The adrenal association is not a brain-cell measurement. The large AlphaGenome DNase prediction remains a **model prediction**, not evidence that the predicted DNA change caused a measured eQTL or phenotype.

The prior Brain6 trait–trait chr3 coloc and the present separate disease–eQTL fits answer different questions; the former does not resolve the latter. None validates a causal LSAMP pathway, repairs the failed LAVA QC, or supplies independent two-trait locus replication. A denser multi-signal analysis with appropriate ancestry-matched LD, fully harmonized source data and experimental validation is needed.

[Protocol](PROTOCOL.md), [QTL receipt](qtl_receipt.json), [matching receipt](match_receipt.json), source-specific matched tables, [conditional H4 variant tables](conditional_h4.tsv), [code](extract_qtl.py) and [verification](verification.json) retain the analysis. The conditional variant probabilities apply only given H4 and should not be read as unconditional causal probabilities, especially where H4 is low. A separate Python implementation reproduced all 48 posterior rows to numerical precision; this is an arithmetic check, not statistical calibration. Indexed GTEx v8 QTLs came from the [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/); remote source and index files were not completely hashed. The original Brain6 regional input hashes and coordinate-liftover chain hash are recorded in the matching receipt.

A post-hoc [reference LD check](lead_ld/RESULTS.md) found that the strongest adrenal LSAMP eQTL marker, rs6804394, has r²≈0.061 with rs62264764 and r²≈0.026 with the matched insomnia lead in 503 European reference people. This supports a follow-up distinction between the observed lead signals but is not conditional analysis or proof of different causal variants.

An [exploratory multi-signal pilot](multisignal_pilot/RESULTS.md) returned one eQTL and one insomnia credible set with low H4, but no ADHD credible set. It is not admitted as validated because study-matched LD and complete source harmonization are unavailable and LD–Z diagnostics show substantial outliers. It does not exclude secondary shared signals.

A [three-trait HyPrColoc pilot](hyprcoloc_pilot/RESULTS.md), assuming independent study statistics, did not return a cluster including adrenal LSAMP expression under the tested settings. It is exploratory and does not rule out secondary sharing.

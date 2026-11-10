# Public latent frailty GWAS: a possible follow-up

[Foote et al.](https://doi.org/10.1038/s41588-025-02269-0) modeled the genetic covariance of 30 frailty-related deficits, creating a general frailty factor and six additional factors. The [GWAS Catalog release](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/) contains genome-wide summary statistics for the general factor (`GCST90624046`), limited social support (`047`), unhealthy lifestyle (`048`), multimorbidity (`049`), metabolic problems (`050`), poorer cognition (`051`) and disability (`052`). The accompanying `053` is pneumonia, not an eighth frailty factor.

These factors are **not** clinical FI or the five-domain Fried score. Their labels are theoretical summaries of latent genetic signals, and the general model includes insomnia as one of its 30 indicators. A factor without insomnia as a named indicator is not equivalent to the full FI with its insomnia item removed: its effects are conditional on the entire fitted genomic model and the general factor. The constituent GWAS use overlapping or potentially overlapping samples, including UK Biobank. Consequently, association with an insomnia dataset does not automatically establish independent replication or quantify the effect of removing a sleep question.

The catalog metadata says the seven latent-factor files use GRCh37 and filtered SNPs with significant factor-specific Q heterogeneity. Their `n` field is the Genomic SEM effective expected sample size **N-hat**, not a participant count. Files are labeled not harmonized. The author-reported N-hat for the general factor is 2,105,711, which must not be described as 2.1 million unique people. The Q-filtered release is a selected subset, so absent variants cannot be interpreted as no association. Before applying LDSC or cross-trait methods, validate whether the factor-specific estimator, N-hat, selection and sample overlap meet that method's assumptions; do not treat these files as ordinary unfiltered single-phenotype GWAS without checking.

## What the residual factors measure

The authors’ supplementary note defines the **multimorbidity factor** through genetic covariance between number of non-cancer illnesses and mean arterial pressure; they interpret it as potentially reflecting cardiovascular-related multimorbidity. The **disability factor** is defined by poor self-rated overall health, long-standing illness/disability/infirmity, and an eye disorder or problem; they interpret it as potentially reflecting chronic disability. These are model-derived genetic factors, not measured clinical diagnoses of those pathways. Insomnia is not among the named indicators defining either residual factor, but both belong to a model whose general factor includes insomnia. Therefore, neither is a sleep-excluded FI, and their associations cannot isolate the effect of removing the insomnia item. See the supplementary note accompanying [Foote et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12626885/).

## Expanded public lookup across all seven factors

I also retrieved indexed results for the four remaining residual factors: limited social support, unhealthy lifestyle, metabolic problems and poorer cognition. Their named defining indicators are respectively partnership/confiding; smoking/arterial stiffness/physical inactivity; body fat/walking pace/basal metabolism/breathlessness; and fluid intelligence/self-rated health/walking pace. These are the authors’ model-derived interpretations, not direct measures of causal pathways. The complete aligned results, standard errors, p-values and Q metrics for all 22 markers and seven factors are in the [seven-factor lookup](foote_seven_factors_shared22_lookup.csv).

| Factor | Exact matched markers | Matched markers with unadjusted p < 0.05 |
|---|---:|---:|
| General frailty | 4 | 3 |
| Limited social support | 20 | 2 |
| Unhealthy lifestyle | 20 | 4 |
| Multimorbidity | 14 | 4 |
| Metabolic problems | 8 | 5 |
| Poorer cognition | 20 | 8 |
| Disability | 20 | 7 |

These exploratory counts depend on different filtered marker sets, correlated variants, post-selection and unadjusted testing. They cannot rank the factors, estimate true sharing, or establish that a candidate acts through a named phenotype. For illustration, rs34290943 near LAMB2 has a nominal association in the poorer-cognition factor (p = 3.23 × 10⁻⁶), while its corresponding disability result is p = 0.572; neither result establishes a mechanism or proves a difference between factors. Formal comparison would require the covariance of factor-effect estimates, careful SNP representation and an appropriately prespecified multiplicity strategy. The four new indexes matched the catalog MD5s, but full compressed files were not downloaded and full-file checksums were not verified. Candidate genetic coordinates and alleles have the same limitations described below.

## Exploratory lookups of the 22 markers

Using the catalog's GRCh38 harmonized, indexed versions, I queried the 22 markers in the general, multimorbidity and disability files. Exact chromosome, position, rsID and allele-pair checks returned **4, 14 and 20** unique rows, respectively. Of the returned rows, 3, 4 and 7 had unadjusted p < 0.05. The files are Q-filtered and have different coverage; these denominators and nominal counts are not a ranking of factors, a replication rate or evidence that missing markers have no effect.

For the five markers previously nominal in all eight inspected GWAS result sets, the available results are:

| Marker; counted allele | General factor p | Multimorbidity factor p | Disability factor p |
|---|---:|---:|---:|
| rs12536335 A | No exact row returned | No exact row returned | 0.0534 |
| rs817509 G | 7.47 × 10⁻⁵ | 0.928 | 0.0379 |
| rs2838788 T | No exact row returned | 0.698 | 0.684 |
| rs34290943 T | No exact row returned | No exact row returned | 0.572 |
| rs12671330 C | No exact row returned | No exact row returned | 0.000110 |

These are post-selected, unadjusted lookups in latent factors, with differing factor definitions and study overlap. They do not demonstrate that variants act on a disability or multimorbidity pathway, nor that a nonsignificant result is a null effect. The two FOXP2-labeled markers are correlated. Insomnia is an indicator of the general model, so the general-factor association is especially unsuitable as independent sleep replication. Raw effects and standard errors on the latent-factor scales are in the [22-marker matrix](foote_three_factors_shared22_lookup.csv); do not add those coefficients to FI, FFS or insomnia coefficients.

The three public tabix indexes matched their catalog MD5s. Exact indexed rows were retrieved remotely, but full-file checksum verification was **not** completed: the original full-file streams timed out, and retries were stopped to avoid appending duplicate data. An empty indexed lookup can reflect selection, harmonization, representation or missingness; the reason is not established. Earlier candidate build conversion is assumed, and exact pair matching does not settle the original multiallelic representation. The per-factor receipts record indexed rows and index hashes; script: `work/foote_factor_indexed_lookup.py`. No participant data were used. Relevant file metadata are in the accession directories, for example [general factor](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/GCST90624046/) and [disability](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/GCST90624052/).

The existing [published multi-score review](multiscore_frailty_prediction_review.md) describes how the authors assessed factor-based scores in FI cohorts. Those published results do not validate our sparse 22-marker count.

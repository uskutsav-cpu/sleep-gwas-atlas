# Cell matched molecular prediction and reserved experimental validation

The response model predicts the magnitude of DNA effects measured in laboratory reporter experiments. The reserved results below determine whether it improves on typical measured effects and whether its uncertainty is useful. Human insomnia and frailty risk are separate outcomes.

## Experimental data and model fitting

The locked development expansion added 41,451 original absolute fitted log2 effect measurements from 13,817 OL27 DNA designs measured in HepG2, K562 and SK-N-SH, giving 654,612 development observations in total. Each connected genomic and shared-DNA component remained in one partition. Original arrays and partitions were preserved as exact prefixes, and all OL33-connected and blood-test regions stayed outside fitting. No detection, significance, activity, effect-size or standard-error cutoff selected examples. Original unsigned targets were retained; unproven coefficient signs and allele activity targets were left missing.

Three fixed author-trained DNA models supplied summaries of each complete source reference/alternate sequence. A small head with 8,461 trainable parameters learned effect magnitude by backpropagation; the author models stayed fixed. One author model per cell was used, rather than the full published ensemble. The selected common epoch was 1, with mean development error improvement of 0.77% over training-only cell/culture medians across seven omitted publications. Whole Siraj omission removed its old and new libraries together. These selected development results do not establish transfer to new experiments.

Final three-seed ensembles used identical observations, shared initial parameters, batches, losses and update counts with and without three added AVI inputs. AVI is the designated permissive composite score, which includes AlphaGenome-related information and other annotations. It was used as an input only. Its contribution cannot be attributed uniquely to AlphaGenome, and absence from every AVI/AlphaGenome training resource has not been established.

## Reserved experimental tests

The independent-publication blood test contained 236 physically verified designs, representing 49 distinct variants in 46 source regions, with measurements in four cell models. Exact K562 was the primary known-line test; HL60, CHRF and THP1 remained unfamiliar-line diagnostics with no authorized predictions. The held OL33 library supplied 39,723 measurements of 13,241 designs in three known lines. It comes from the existing Siraj publication and lab, so it is not independent-publication or donor replication. A further 1,932 OL27 reserved measurements are descriptive.

Full actual 148 bp blood and 200 bp library sequences were used without padding or invented genomic context. Canonical genomic neighborhoods and literal shared DNA were screened against the published libraries covering all actual author pretraining sequences. Unmapped source coordinates limit wider-neighborhood claims. Published paired GLM coefficients and source count aggregates were checked, but the original raw barcode pipeline and GLM were not independently rebuilt.

Both ensembles, uncertainty models, acceptance rules and predictions were frozen before the numeric target join. Accuracy used absolute original tile-level Log2Skew. The primary error and correlation measures give equal weight to cell lines; correlations are computed within each line before averaging. Intervals below are empirical targets, with 2,000 paired bootstrap samples of fixed one-megabase windows for accuracy comparisons. These windows are not ancestry-specific LD blocks.

| Test | Input | Mean absolute error | Typical effect error | Relative error improvement | Within line correlation | Interval coverage | Accepted before outcomes |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Independent blood K562 | sequence | 0.140531 | 0.136613 | -2.87% | -0.0872 | 99.58% | 0 / 236 |
| Independent blood K562 | AVI | 0.141287 | 0.136613 | -3.42% | -0.0739 | 99.58% | 0 / 236 |
| Held OL33 library | sequence | 0.171032 | 0.175737 | 2.68% | 0.0673 | 94.86% | 0 / 39723 |
| Held OL33 library | AVI | 0.170888 | 0.175737 | 2.76% | 0.0674 | 94.86% | 0 / 39723 |

Negative improvement means worse prediction than the typical-effect baseline. Complete confidence intervals, per-line results and direct nearby/cell difference comparisons are retained in reserved_outcome_evaluation.json. The independent audit separately reconstructed all original targets and primary error, correlation, coverage and bootstrap results; pair point estimates were independently checked, but pair confidence intervals were not independently replayed.

## Confidence and nearby cell differences

Development checks excluded each outer publication from base fitting, confidence fitting and interval scaling. Confidence training used other publications original calibration regions from that same omitted-source base model. The original source-protected Siraj check still achieved only 76.43% overall coverage, with all cases withheld. Final scaling used conservative source-specific residual quantiles and absolute error floors. Exact unfamiliar lines, unusual inputs and cases without sufficient predicted benefit or narrow uncertainty were withheld. High coverage achieved by broad intervals, or withholding every case, is not useful prediction.

sequence: the held OL33 test provided 1,119 nearby comparisons and 39,723 comparisons of the same physical design across cell lines. Direct difference error, median-difference baselines, within-fixed-cell-pair correlations and interval coverage are reported separately. A pooled difference correlation alone does not establish variant-specific cell effects.

AVI: the held OL33 test provided 1,119 nearby comparisons and 39,723 comparisons of the same physical design across cell lines. Direct difference error, median-difference baselines, within-fixed-cell-pair correlations and interval coverage are reported separately. A pooled difference correlation alone does not establish variant-specific cell effects.

## AlphaGenome related input decision

Retain AVI under the prespecified tests: **False**. Both primary relative-gain and confidence-interval gates passed: **False**. Both absolute accuracy gates passed: **False**. Local OL33 confirmation was sufficiently large: **True**, and ranking preservation passed: **False**.

## Manuscript readiness

Completed evidence comprises source and allele/physical-DNA audits, reproducible training, matched added-input testing, source-protected confidence checks and reserved experimental evaluation with an independent numerical audit. Development gains remain exploratory; original V7 and V9 adverse tests retain their original status. The experimental model supplies no independent person-level genetic-score validation, sleep/fatigue-excluded frailty GWAS, genetic correlation, disease causation or causal variant-to-gene proof. The Nature Communications genetic manuscript remains unready while those scientific validation gaps remain.

The reserved outcomes are now revealed and cannot be reused as untouched tests of a subsequently changed model. Any further model development must preserve these findings and obtain new appropriate validation before claiming improved transfer.

Primary sources: [Siraj et al.](https://doi.org/10.1038/s41586-026-10121-6), [Tardaguila et al.](https://doi.org/10.1016/j.celrep.2026.116937), and the pinned author pretrained DNA model/library source records retained in next_source_screen.

### Nearby changes and differences between cell lines

In the held library, the sequence model correctly ordered absolute measured effects in 52.73% of 1,119 nearby comparisons across 181 one-megabase windows. The bootstrap interval for ordering was 49.38% to 55.82%, including chance ordering. Its paired effect-difference error improvement over predicting no nearby difference was 0.51%, with an interval crossing zero. The independent blood publication supplied only two eligible nearby comparisons, which cannot establish local replication.

For 39,723 comparisons of the same physical designs between the three cell lines, ordering accuracy was 43.30%. The average correlation within fixed pairs of cell lines was -0.0153. Predicting no between-cell difference gave lower absolute error (0.175347) than the model (0.192371). Improvement over the less accurate training-median difference therefore does not establish variant-specific cell learning. These results do not support cell-specific prioritization of sleep or frailty variants.

Adding AVI increased independent blood prediction error by 0.54% relative to sequence alone; the paired window bootstrap interval for the error increase was 0.03% to 1.12%. The small held-library gain did not satisfy the predefined retention rule. The composite input remains excluded from authorized biological predictions.

The nominal 90% intervals achieved high coverage with a median halfwidth of about 0.45 log2 units, while both arms withheld all cases. This is conservative withholding without useful answer availability. The accuracy requirements were a minimum 5% error improvement over both training-only typical-effect baselines, positive lower confidence bounds and a within-cell-line correlation of at least 0.1 with a positive lower bound. The added-input requirements also included at least 1% improvement in both primary tests and preservation of nearby ranking.

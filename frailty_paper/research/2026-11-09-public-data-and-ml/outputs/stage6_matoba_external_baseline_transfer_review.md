# Fixed sequence baseline on a second neuronal reporter study

The existing Retallick-trained 1–3-mer logistic baseline was reconstructed using its original fixed training and scaling code, then evaluated without fitting or tuning on Matoba et al.'s vehicle-condition allelic-effect results. We restricted Matoba records to chromosomes 1 and 2, excluded from the original model's training, and to the 108 unique sequence pairs passing the GRCh37 reference-position and allele-pair check in the [source audit](stage6_matoba_neural_reporter_benchmark_feasibility.md). This is a post hoc cross-study methods test, not a participant-level prediction or validation of the project's genetic scores.

| Matoba sequence orientation | AUROC | Average precision |
|---|---:|---:|
| As released | 0.533 | 0.630 |
| Reverse-complemented to reference-plus | 0.547 | 0.659 |

The 108 reported allelic tests comprise 64 adjusted-p < 0.05 calls and 44 nonsignificant calls; detected-call prevalence is 0.593. Constant-score AUROC is 0.5 and average precision is the prevalence. The two orientations were assessed on the same variants as a sensitivity, not independent replications; the small apparent difference is not a formal statistical comparison. Both AUROCs are near chance and these results do not establish useful generalization. The precision–recall baseline differs from the earlier Retallick test because the selected Matoba subset has a much higher positive-call prevalence.

The target is a study-specific statistical detection call, **not** a ground-truth functional/nonfunctional classification; nonsignificance can reflect low power. Matoba's allelic-effect table is a selected subset of its activity-tested constructs, and the two studies differ in cell models, sequence lengths (150 vs 199 bases), designs, donors and underlying variant selection. Entire chromosomes were held out from baseline training. An exact-sequence check found no full-length matches, but homologous sequence overlap, shared donors and other source/pretraining overlap have not been fully assessed. The same test was chosen after inspecting its availability, and 108 results do not establish reproducible performance. No AlphaGenome outputs or derivatives were used. The paper describes CC BY-NC-ND 4.0; rights for broader dataset/model reuse require separate verification, and only aggregate metrics are reported here.

Reproduction: bundled Python running `work/matoba_fixed_baseline_transfer.py`, which rebuilds the same fixed prior model, uses existing source workbooks and `stage6_matoba_reference_alignment.json`, evaluates both orientations and records [aggregate metrics](stage6_matoba_external_baseline_transfer.json). [Original model and uncertainty](stage6_neuronal_mpra_ml_baseline_review.md); [source and reference-alignment audit](stage6_matoba_neural_reporter_benchmark_feasibility.md). Clinical prediction and independent score validation remain outstanding.

## Uncertainty accounting for reported loci

The 108 evaluated records map to 26 author-reported locus labels; the largest locus contains 15 variants. Keeping the model and reference-plus orientation fixed, a 5,000-draw bootstrap that resamples these locus labels with replacement yielded AUROC **0.547**, with a 2.5th–97.5th percentile range of **0.430–0.661**. This illustrates substantial uncertainty and is compatible with chance-level ranking; it does not prove equivalence to chance or rule out useful performance in other settings.

The interval treats the 26 reported loci as exchangeable independent sampling units, an approximation that may not reflect true biological or genomic dependence; loci were selected for the original study and only chromosomes 1–2 were assessed. It includes neither model-training uncertainty nor assay, label, selection, donor or license uncertainty. The bootstrap samples whole reported loci without class stratification, preserving within-locus records but allowing class prevalence to vary. All 5,000 draws included both classes. The label-based partition is not a guarantee of independent genomic blocks. No individual-level clinical data or AlphaGenome predictions were used.

Reproduce with bundled Python and `work/matoba_locus_bootstrap.py`; [aggregate receipt](stage6_matoba_locus_bootstrap.json) records the fixed seed and draw counts. The point AUROC is independently recomputed with a tie-aware pairwise formula and checked against the original transfer output.

## Exact construct-sequence overlap check

An additional post hoc check compared both 150-base alleles of each of the 108 Matoba test pairs, in both orientations, against both 199-base constructs of all 801 Retallick training pairs. **Zero test pairs had an exact full-length substring match.** This addresses exact reuse of the tested constructs only; shorter or approximate matches, homologous regions, donor overlap and other indirect data overlap remain untested. It does not establish statistical independence or validate the genetic scores. Reproduce with `work/matoba_training_sequence_overlap.py`; the [aggregate receipt](stage6_matoba_training_sequence_overlap.json) contains no individual sequences.

## Orientation sensitivity of variant rankings

The similar AUROCs for the two orientations conceal substantial changes to the individual rankings. Across the same 108 pairs, the Spearman correlation between tie-averaged ranks was **−0.113** and the mean absolute rank shift was **37.9 positions**. The median absolute difference in model output was 0.064; these outputs are not calibrated probabilities. These are descriptive, post hoc results on a selected reporter dataset, not a clinical score. They show that this particular k-mer baseline is sensitive to sequence orientation and is not ready for stable variant prioritization without a prespecified orientation policy or a tested reverse-complement-invariant design.

For the paired AUROC difference (reference-plus minus as-released), the point estimate was +0.014 and a 5,000-draw reported-locus bootstrap percentile interval was **−0.169 to +0.216**. This does not establish that either orientation performs better. The interval assumes the reported locus labels are exchangeable independent units, keeps the trained model fixed, and does not capture training, source-selection or donor uncertainty. No new model was fitted and no variant-level predictions are distributed. Reproduce with `work/matoba_orientation_sensitivity.py`; [aggregate result](stage6_matoba_orientation_sensitivity.json).

## Exploratory strand-invariant baseline

As a direct response to the orientation sensitivity, we fitted a new logistic baseline using the original training chromosomes only. For each allele pair, we average absolute 1–3-mer count-difference features with their reverse-complement features before training-only scaling. The model uses the original learning rate, penalty and iteration count. On all 108 Matoba test pairs, the features were exactly identical when both allele sequences were reverse-complemented. This establishes invariance for this implementation, not predictive validity.

| Dataset and held-out chromosomes | Calls | Strand-invariant AUROC | Average precision |
|---|---:|---:|---:|
| Retallick 1–2 | 161; 32 detected | 0.482 | 0.247 |
| Matoba 1–2 | 108; 64 detected | 0.569 | 0.675 |

The original orientation-sensitive baseline achieved Retallick AUROC 0.517 and Matoba AUROC 0.533–0.547 depending on orientation. This redesign was made **after inspecting both test sets and the orientation problem**, so these numbers are exploratory comparisons, not evidence of improvement or independent confirmation. The studies use different assays and selected variant sets; average precision also reflects different positive-call prevalence. The target remains a statistical detection call, not functional truth. No untouched dataset, patient outcomes or AlphaGenome predictions were used. Freeze the transformation and evaluation before any further assessment on an untouched suitably licensed independent dataset; neither model is ready for score validation or clinical prediction.

Reproduce with `work/matoba_invariant_baseline.py`; the [aggregate receipt](stage6_matoba_invariant_baseline.json) includes an exact feature-invariance assertion. Source-data licensing caveats above apply; no learned weights or sequence-level outputs are distributed.

### Paired uncertainty for the redesign

On the same Matoba pairs, the invariant-minus-original reference-plus AUROC difference is **0.023**. A paired 5,000-draw bootstrap by the 26 reported loci gives a percentile interval of **−0.079 to +0.103**, consistent with both a decrease and an increase. This interval is conditional on the fixed models and assumes exchangeable, independent reported loci; it excludes model-training and dataset-selection uncertainty. Because the change was designed after examining this test set, the comparison cannot establish improvement even if the interval were narrower. Reproduce with `work/matoba_invariant_paired_bootstrap.py`; [aggregate results](stage6_matoba_invariant_paired_bootstrap.json).

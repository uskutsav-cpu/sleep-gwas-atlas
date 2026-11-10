# AlphaGenome-assisted model development, version 2

This plan precedes v2 feature fitting and selection. Preserve v1 results. The goal is the strongest defensible research model supported by the available free data, not a guaranteed positive result or a clinical predictor.

## Allowed AlphaGenome input

Use only the downloadable AVI SNV raw/PHRED scores explicitly designated under the Permissive Use Downloadable Artifacts section on the official downloads page, observed 6 November 2026. The current service terms' training restriction explicitly excludes those artifacts. The detailed live API predictions and AVI feature breakdowns are excluded. The exact-allele indexed mirror is sourced from the Genomic Resource Repository; verify the index MD5, retain resource metadata, exact coordinates and alleles, retrieval status and hashes. AVI is a composite predicted molecular-impact score (including conservation and other annotations), not measured health or laboratory truth. Do not call this learning directly from AlphaGenome weights.

Official sources: https://deepmind.google.com/science/alphagenome/terms and https://deepmind.google.com/science/alphagenome/downloads . The pages were inspected in the signed-in research browser without accepting new terms. The downloaded-AVI exception does not authorize training from all other outputs. No clinical use.

## Measurements and honest evaluation

The primary supervised target remains the measured Salomon NGN2 allelic reporter call. Use the v1 audited 37,349 unique autosomal sequence pairs and its chromosome splits: 1–16 training, 17–18 validation, 19–22 reused comparison benchmark. The v1 final-test results have already been inspected: this last split is NOT a new independent test and must be labelled as a reused benchmark. Do not tune on its new results. Model selection uses validation average precision only; preserve every attempted configuration. Seek an additional separately collected assay for evaluation, checking sequence/region overlaps, genome builds and assay differences. Such testing is not independent genetic-score validation.

Inputs: paired DNA features, the pinned original HyenaDNA frozen representations, and permitted AVI scores. Missing AVI scores remain missing with an availability indicator; never substitute a biological zero. The detailed API outputs, GWAS p-values, gene/disease labels, source outcome-derived features and reporter reference/alternate activities are excluded from inputs. Only measured target labels enter supervised loss. AVI's own model training/benchmark overlap is unresolved and must be disclosed. Prior exposure to the public human genome is not excluded by the supervised split.

## Fixed first search

1. AVI score alone, both its direct ranking and a supervised recalibration.
2. Regularized logistic models with DNA counts alone versus the same counts plus AVI (C 0.01, 0.1, 1, 10).
3. Histogram gradient boosting with and without AVI: depths 2, 3, 5; learning rates 0.03, 0.08; 200 trees, minimum leaf 40, L2=10. Use explicit chromosome validation instead of random internal early stopping.
4. Small neural heads on original frozen HyenaDNA features, with and without AVI. Hidden sizes 32 and 128, three seeds, AdamW, at most 50 epochs, validation-AP stopping after 7 nonimproving epochs. Use paired settings to distinguish input improvement from architecture changes. This is real backpropagation through the head; v1 already performed encoder fine-tuning.
5. If AVI improves validation, test a sequence student initialized from v1's pretrained encoder with an auxiliary AVI regression target. Use only training chromosome AVI targets; measured reporter labels remain the primary objective. This separate experiment tests imitation as an auxiliary task, not independent experimental truth. Log all amendments before accessing benchmark results.

All preprocessing fits on training chromosomes only. Select the overall recommendation by validation AP, including the v1 simple comparator. Report ROC AUC, AP, prevalence, precision/recall and calibration qualifications. Model outputs are uncalibrated research scores. A model may be strong in one assay and weak elsewhere. Do not report ensemble or candidate performance without identifying the selected model and its validation rule. Candidate scores do not establish disease associations or gene/cell mechanisms.

## Required artifacts

Exact-allele AVI retrieval audit; all search configuration/validation results; saved selected weights; comparable ablations; inference/reload checks; training curves; repeated-run variability for neural heads; results with clear separation of validation, reused benchmark and external assay; a research model card and reproducible code. Additional validation and an evidence-based stopping decision remain required before calling the full goal complete.

## Auxiliary-student configuration locked before v2 validation results

The conditional native-weight student uses seed 1207, the original pretrained encoder and the v1 frozen-head initialization. Its assay head sees DNA features only. A separate linear head predicts the standardized/clipped raw AVI target; the joint loss is weighted measured-call BCE plus 0.05 times masked Huber teacher loss. Missing teacher values contribute no teacher loss. Warm up heads for at most 20 epochs (patience 5), then update the last native block and final normalization for at most 4 epochs (patience 2), batch 32, encoder LR 1e-5 and head LR 3e-4. Choose across both phases by measured-assay validation AP. Native updates that fail to improve the frozen student must be reported, not silently promoted. This is additional learning from permitted AVI, not proof that an AVI prediction is a true laboratory effect. The exact code is train_student.py. No reused benchmark or external model performance was available when this configuration was written.

## Separate-study evaluation amendment, before any v2 model performance

Retallick-Townsley baseline constructs were independently reconstructed from Supplementary Data 2 and counts cross-checked against the existing project source audit. All 962 199-base constructs match hg38, and 30 lie within 1 kb of any Salomon contrast. Reserve the remaining 932 pairs (134 detected calls) for evaluation only. Input native 270-base contexts are checked against the 199-base constructs; extra flanks were not part of the physical reporter assay. Source labels are not used to select these exclusions. The study's broader project aggregates and an earlier baseline were previously inspected, and model pretraining/AVI exposure is unresolved; describe the check as external to current fitting rather than an entirely unseen or foundation-model-independent study. Its sequence derivatives remain local because the article states CC BY-NC-ND 4.0.

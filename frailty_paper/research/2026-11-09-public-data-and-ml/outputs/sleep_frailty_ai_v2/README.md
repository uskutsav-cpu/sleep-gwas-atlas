# AlphaGenome-assisted Sleep–Frailty research models, v2

The selected model is **logistic_kmer_C10**, chosen by the largest average precision (AP) on Salomon chromosomes 17–18. All 43 tested model configurations are retained. A larger model or an AlphaGenome input is not assumed to be better. The separate assay and reused benchmark never enter the selection rule. The finite multi-model search used only 78 positive validation calls; validation optimism is an important limitation.

The auxiliary native-weight experiment was extended after observing round1 benchmark and separate-study results, despite its original validation trigger not firing. The architecture/loss settings were already specified, but this execution decision is exploratory. Round1 is preserved under round1/. Both existing evaluation sets are reused development comparisons for this extension; a further separate source is required for a fresh evaluation.

## What AlphaGenome contributes

The training input is the downloadable AVI SNV raw/PHRED score, not the detailed API predictions. The official [downloads page](https://deepmind.google.com/science/alphagenome/downloads) lists AVI SNV scores as a permissive-use artifact. The [service terms](https://deepmind.google.com/science/alphagenome/terms), inspected on 6 November 2026, explicitly except those artifacts from the similar-model-training restriction. This does not authorize training from the other API outputs or AVI feature breakdowns. Exact-allele indexed retrieval from the [Genomic Resource Repository](https://grr.iossifovlab.com/hg38/scores/AVI/index.html) found 37,382 scores. Missing scores remain unavailable, not zero biological effects. The mirror index checksum was verified; the whole 88 GB table was not independently hashed.

AVI combines predicted functional information and other annotations; it is neither a measured neuronal effect nor a sleep/frailty label. Its model-training and benchmark exposure to the evaluated assays has not been established. Adding AVI can transfer information from a larger model, but does not create new experimental validation.

The paired search compares identical fitting settings with and without AVI: regularized logistic models, tree boosting and small neural heads on frozen pretrained DNA representations. Neural heads were optimized by backpropagation across three seeds. A separate native-weight student was conditionally tested if the AVI-assisted candidates improved validation. Student status: **trained**. Its primary target is the measured reporter call, and AVI is an auxiliary target without being an input to its assay prediction head. Learning an AVI target is imitation; only measured-assay evaluation can assess experimental usefulness. See `student_config.json` and `student_status.json` when trained.

## Data and evaluation

Primary fitting uses the v1 audited 37,349 Salomon human excitatory-neuron 270-base DNA pairs: 31,333 on training chromosomes 1–16, 2,474 validation pairs on 17–18, and 3,542 pairs on chromosomes 19–22. The last group was already inspected in v1. It is a **reused comparison benchmark**, not a fresh independent test. V1 results and weights are preserved. Only paired DNA and permitted AVI enter new model inputs; gene labels, disease GWAS statistics, detailed AlphaGenome outputs and outcome-derived source annotations do not.

For a separate-study check we reconstructed Retallick-Townsley et al.'s baseline neuronal reporter assay ([paper](https://doi.org/10.1038/s41467-026-77114-x)). All 962 source DNA constructs matched hg38 exactly, including one observed allele. Removing variants within 1 kb of any Salomon variant retained 932 pairs with 134 detected calls. Public native 270-base contexts were used for consistency with model input length, after confirming the measured 199-base constructs inside them. The extra native flanks were not tested by that experiment, and biological state/assay selection differ. Region exclusion does not prove no sequence homology or foundation-model exposure. This study was previously inspected elsewhere in the project; it is external to current fitting, not an entirely unseen publication. No external labels train or select these models.

## Results

| Model | Validation AP | Reused benchmark AUC | Reused benchmark AP | Separate assay AUC | Separate assay AP |
|---|---:|---:|---:|---:|---:|
| logistic_kmer_C10 **selected** | 0.0536 | 0.520 | 0.0465 | 0.536 | 0.1590 |
| logistic_kmer_C1 | 0.0532 | 0.524 | 0.0499 | 0.534 | 0.1568 |
| v1_simple_baseline | 0.0498 | 0.583 | 0.0487 | 0.481 | 0.1496 |
| logistic_kmer_avi_C10 | 0.0497 | 0.559 | 0.0455 | 0.545 | 0.1575 |
| logistic_kmer_avi_C1 | 0.0489 | 0.562 | 0.0464 | 0.544 | 0.1570 |
| logistic_kmer_C0.1 | 0.0478 | 0.532 | 0.0564 | 0.530 | 0.1538 |
| logistic_kmer_avi_C0.1 | 0.0477 | 0.572 | 0.0522 | 0.538 | 0.1551 |
| logistic_kmer_C0.01 | 0.0470 | 0.543 | 0.0523 | 0.524 | 0.1548 |
| logistic_kmer_avi_C0.01 | 0.0467 | 0.577 | 0.0494 | 0.536 | 0.1556 |
| neural_hyena_avi_w128_s101 | 0.0426 | 0.600 | 0.0412 | 0.478 | 0.1515 |
| ensemble_hyena_avi_w128 | 0.0421 | 0.597 | 0.0400 | 0.487 | 0.1533 |
| neural_hyena_avi_w128_s73 | 0.0412 | 0.593 | 0.0394 | 0.478 | 0.1482 |
| neural_hyena_avi_w128_s42 | 0.0408 | 0.581 | 0.0387 | 0.506 | 0.1581 |
| neural_hyena_avi_w32_s101 | 0.0406 | 0.585 | 0.0398 | 0.515 | 0.1502 |
| neural_hyena_w128_s101 | 0.0404 | 0.569 | 0.0357 | 0.490 | 0.1438 |
| neural_hyena_w128_s73 | 0.0397 | 0.563 | 0.0342 | 0.490 | 0.1443 |
| neural_hyena_w32_s101 | 0.0396 | 0.542 | 0.0324 | 0.491 | 0.1551 |
| ensemble_hyena_w128 | 0.0394 | 0.565 | 0.0362 | 0.491 | 0.1492 |
| neural_hyena_w128_s42 | 0.0392 | 0.542 | 0.0358 | 0.510 | 0.1478 |
| ensemble_hyena_avi_w32 | 0.0388 | 0.590 | 0.0423 | 0.502 | 0.1474 |
| neural_hyena_w32_s73 | 0.0382 | 0.555 | 0.0350 | 0.491 | 0.1372 |
| neural_hyena_avi_w32_s42 | 0.0375 | 0.594 | 0.0409 | 0.496 | 0.1499 |
| boost_kmer_avi_d5_lr0.03 | 0.0373 | 0.518 | 0.0449 | 0.511 | 0.1453 |
| logistic_avi_C10 | 0.0372 | 0.586 | 0.0426 | 0.540 | 0.1583 |
| avi_direct_rank | 0.0372 | 0.586 | 0.0426 | 0.542 | 0.1584 |
| logistic_avi_C0.01 | 0.0372 | 0.586 | 0.0426 | 0.542 | 0.1584 |
| logistic_avi_C0.1 | 0.0372 | 0.586 | 0.0426 | 0.542 | 0.1584 |
| logistic_avi_C1 | 0.0372 | 0.586 | 0.0426 | 0.541 | 0.1584 |
| neural_hyena_w32_s42 | 0.0370 | 0.561 | 0.0348 | 0.495 | 0.1568 |
| boost_kmer_d2_lr0.03 | 0.0369 | 0.504 | 0.0432 | 0.489 | 0.1478 |
| neural_hyena_avi_w32_s73 | 0.0366 | 0.578 | 0.0454 | 0.500 | 0.1420 |
| ensemble_hyena_w32 | 0.0365 | 0.554 | 0.0346 | 0.490 | 0.1525 |
| boost_kmer_avi_d2_lr0.03 | 0.0364 | 0.536 | 0.0535 | 0.514 | 0.1433 |
| avi_auxiliary_student | 0.0359 | 0.529 | 0.0333 | 0.515 | 0.1500 |
| boost_kmer_avi_d3_lr0.08 | 0.0357 | 0.534 | 0.0509 | 0.496 | 0.1370 |
| boost_kmer_d3_lr0.03 | 0.0356 | 0.493 | 0.0359 | 0.487 | 0.1477 |
| boost_kmer_d2_lr0.08 | 0.0356 | 0.518 | 0.0459 | 0.483 | 0.1400 |
| boost_kmer_avi_d3_lr0.03 | 0.0354 | 0.534 | 0.0511 | 0.507 | 0.1414 |
| boost_kmer_avi_d2_lr0.08 | 0.0347 | 0.549 | 0.0578 | 0.511 | 0.1416 |
| boost_kmer_d5_lr0.03 | 0.0347 | 0.487 | 0.0298 | 0.502 | 0.1493 |
| boost_kmer_d3_lr0.08 | 0.0339 | 0.507 | 0.0436 | 0.489 | 0.1469 |
| boost_kmer_avi_d5_lr0.08 | 0.0339 | 0.536 | 0.0464 | 0.494 | 0.1436 |
| boost_kmer_d5_lr0.08 | 0.0323 | 0.505 | 0.0365 | 0.478 | 0.1429 |

AUC 0.5 denotes chance ranking, not 50% clinical accuracy. AP depends on positive prevalence; compare models within the same assay rather than comparing AP values across different datasets. Binary calls reflect statistical detection versus nonsignificance, not proven biological active/inactive truth. Class weighting and source selection mean raw neural outputs are uncalibrated. Direct AVI is ranked through a monotonic transform solely to keep score ranges finite; its Brier statistic is not a probability-calibration result.

`conditional_uncertainty.json` reports paired chromosome-resampling percentile ranges for the selected model against the DNA-only model chosen by validation and against direct AVI. These are conditional descriptive sensitivities, not complete uncertainty guarantees. Only four reused benchmark chromosomes exist. Training, search/selection, study sampling and foundation-model exposure uncertainty are not included. Point improvements alone are insufficient to establish that adding AlphaGenome reliably helps.

The selected model's reused benchmark AUC is 0.520, AP 0.0465. Its separate-study AUC is 0.536 and AP 0.1590. See all model results, including failures to improve. No recommendation was changed to whichever model won the separate-study comparison.

## Verification, scope and remaining work

Frozen encoder caches were replayed on sampled rows against the pinned original weights with exact agreement; source/checkpoint/feature hashes are retained. All selected heads and comparators are reloaded and checked against recorded validation predictions in `checkpoint_replay.json`. Independent R ranking calculations are compared with the Python metrics in `independent_metric_agreement.json`. This research code uses the previously checked custom NumPy/MLX HyenaDNA port; official PyTorch-executable parity remains untested.

The 33 project candidate predictions are exploratory model outputs without matching experimental candidate labels. They do not establish disease association, causality, a causal gene, cell specificity across tissues, genetic correlation, or validation of 15-/22-variant scores. None of these models predicts an individual's insomnia, frailty or circadian disorder. Stronger claims need additional independent assay data, appropriate uncertainty and cell-state testing, and separate genotype/outcome evidence for clinical prediction.

Source sequence-level extracts and external predictions remain local. The Retallick-Townsley article is CC BY-NC-ND 4.0; do not redistribute adapted sequence-level material. Salomon source/data and HyenaDNA licenses remain as documented in the v1 report. Preserve upstream notices. No model package has been published externally.

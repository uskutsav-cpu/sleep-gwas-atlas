# Sleep–Frailty AI v1: a trained neuronal regulatory-effect model

## What was built

We adapted the public **HyenaDNA tiny** pretrained DNA model, which has approximately 450,000 stored scalars, using real measurements from a human excitatory-neuron reporter experiment. The trained model takes a reference DNA sequence and an alternative sequence and returns a score for whether the experiment would detect a difference in their regulatory activity. The target is the study's detected-versus-nonsignificant call, which depends on assay power as well as biology. A nonsignificant call does not prove an inactive variant.

This is a **molecular research prediction model**. It does not predict whether a person will develop insomnia, frailty, or a circadian disorder, and it does not validate the project's 15-/22-variant scores. It cannot identify a causal gene simply from a sequence score. Endogenous regulation, other tissues, aging states and human outcomes require separate validation.

## Actual training

The original pretrained sequence encoder has two blocks, 128 features per DNA position and a 1,026-token capacity. We kept its full 270-base assay contexts, paired reference/alternate sequences, and averaged both DNA strands. Pair features combine sequence context and absolute differences in sequence representations. A 32-unit nonlinear prediction head was trained. Then the last pretrained block and final layer normalization were updated through **backpropagation**, using a smaller learning rate than the head. The first block and nucleotide embedding were frozen.

Native pretrained parameters selected for updates: **217,152**. New head parameters: **16,449**. **25/25 selected native parameter arrays changed** in the saved fine-tuned checkpoint, with their individual changes retained in `test_results.json`. Frozen-head checkpoint selection chose epoch 1; fine-tuning checkpoint selection chose epoch 3. We saved both checkpoints and the simple baseline. A fixed seed, learning rates, class-weighted loss and validation-only selection rules were recorded before final testing.

The run used MLX 0.22.1 on the local Apple GPU. Disk pressure initially prevented installing the runtime; 74 temporary byte-identical public NHANES source copies were replaced with hash-verified links to existing originals, preserving their contents and paths. A feature-extraction attempt was interrupted before any checkpoint or test result existed, and the training was restarted after confirming no live process remained. The final run uses explicit GPU memory/cache limits and reusable frozen-feature caches. No model was redesigned in response to the final test results.

Because the frozen encoder did not improve on the simple comparator in validation, a fixed small convolution model was added before final-test inspection as a stronger local-sequence comparator. It has 23,617 parameters, three shared DNA convolution layers and a small head; all its weights were trained from scratch through backpropagation. It is not pretrained. An initial attempt saturated its prediction head and gave constant validation scores. That attempt was stopped, retained separately, and restarted with parameter-free per-example RMS feature normalization; no cross-example test statistics enter that normalization. The corrected model's selected checkpoint was epoch 2. This addition and debugging reused the validation set and were documented in `CNN_COMPARATOR_PLAN.md`; they were not part of the initial fixed comparison list. All model-fitting and selection choices were completed without inspecting final-test performance.

## Data quality and holdout

The source table had 38,968 rows. After removing duplicated genomic contrasts and nonautosomal/unusable entries, **37,349 unique autosomal pairs** remained. Each pair has 270 valid DNA bases, one verified substitution, a documented coordinate/strand, and a measured readout with an unambiguous source call.

- Training: **31,333** pairs, 809 detected calls, chromosomes 1–16.
- Validation for choosing checkpoints: **2,474** pairs, 78 detected calls, chromosomes 17–18.
- Final supervised test: **3,542** pairs, 100 detected calls, chromosomes 19–22.

Whole chromosomes were separated instead of randomly mixing nearby variants between training and test. Exact reference/alternate sequence matches, including reverse complements, were checked across splits. This is not a guarantee against homologous sequences or study-specific bias. The reference genome, including test chromosomes, may have been seen during unsupervised pretraining; only supervised assay labels are held out. The assay is one selected study, not an independent disease cohort, and earlier project baselines used other parts of this resource. Selection of source variants partly used an existing sequence predictor.

No AlphaGenome output, Enformer prediction, GWAS result, associated-gene label, or outcome-derived annotation entered the new model's features. The simple comparator uses only fixed DNA-pattern counts. Source annotation columns were not inadvertently treated as independent training labels.

## Final test results

| Model | Test ranking AUC | Test average precision |
|---|---:|---:|
| Constant score | 0.500 | 0.0282 |
| Simple DNA-pattern model | 0.583 | 0.0487 |
| Pretrained DNA model + trained head | 0.515 | 0.0296 |
| Fine-tuned DNA model | 0.532 | 0.0339 |
| Small local DNA convolution model | 0.490 | 0.0288 |

The test detected-call prevalence was **2.82%**, which is the random-ranking average-precision reference. A ranking AUC of 0.5 corresponds to chance pair ordering; it is not percent clinical accuracy. Average precision describes enrichment of detected calls toward the top of the ranking and depends strongly on prevalence. It must not be read as a health-risk probability.

The **validation-only** model recommendation is **Simple DNA-pattern model**. Its final-test AUC was 0.583, and average precision was 0.0487. This recommendation did not choose the largest test result. Its validation-selected threshold yielded test precision 4.3% and recall 35.0%; these values concern the selected reporter experiment, not clinical benefit. Raw classifier outputs are uncalibrated because of class weighting and assay selection.

`test_uncertainty.json` gives paired chromosome-bootstrap percentile ranges and comparisons. Only four test chromosomes exist, so these ranges are a rough dependence sensitivity; they exclude training uncertainty, model/split choices, assay differences and the selected tested universe. They do not establish robust performance in a new study or patient population. All five model results, including failures to improve, are retained.

## Verification and implementation qualifications

The compact encoder is a reviewed NumPy/MLX implementation of the pinned public model formulas. Before fitting, its outputs were compared with a separate NumPy/Autograd implementation, long convolutions were checked against direct convolution, and gradients through the encoder were checked against finite differences and Autograd. Allele-swap and reverse-complement symmetry checks passed. The installed runtime's real-FFT derivative failed an identity check, so training instead used complex FFT with a real inverse output; its identity and convolution derivatives passed. See `port_check.json` for actual errors.

Official PyTorch-executable parity has **not** been tested. The checks support this research implementation but do not make it an official LongSafari release. Model checkpoint reload, CPU/GPU prediction agreement and rejection of malformed inputs are checked separately in `inference_check.json` when present.

## Applying it to the project's variants

Public UCSC hg38 reference contexts were obtained for all 33 distinct project candidates. Each context was checked for length, coordinate and reference-base membership in the mapped candidate allele pair. These are public reference sequences, not participant genotypes. `project_reference_context_receipt.json` retains exact sources and hashes. Predictions, when present in `project_candidate_predictions.csv`, are **exploratory reporter-effect scores**, not association evidence or risk estimates. None of these 33 exact contrasts was present in the source assay's earlier candidate screen, so candidate-level accuracy is unknown. The source disease-study multiallelic contrasts and original FFS conversion remain unresolved by this model.

## Sources, permissions and attribution

- Base weights and model card: [LongSafari HyenaDNA tiny](https://huggingface.co/LongSafari/hyenadna-tiny-1k-seqlen-hf), pinned revision `e8c1effa8673814e257e627d2e1eda9ea5a373f6`; model card states BSD-3-Clause. The upstream HyenaDNA implementation repository separately publishes an Apache-2.0 license, retained in the package notices. The custom port and new head are research adaptations; preserve upstream notices on reuse.
- Measured labels and DNA pairs: Salomon et al., *Massively parallel characterization and predictive modelling of neuronal regulatory variation*, [preprint DOI](https://doi.org/10.64898/2026.07.16.738760), [source table and analysis repository](https://github.com/kircherlab/80K-Analysis/tree/a686f1d552f97ebcc0224057b5a835002f7a3dc8). Source software is MIT. The article XML obtained from Europe PMC states CC BY 4.0; article licensing does not automatically resolve rights for every third-party annotation in the repository. This run uses only the authors' DNA pairs and measured calls. The training table and third-party model annotation columns are not redistributed in the model package.
- Reverse-mode numerical check: [HIPS Autograd](https://github.com/HIPS/autograd), version 1.8.0, MIT.

The files are for research review. A stronger biological claim needs an independently collected assay, an endogenous variant–gene test in a relevant cell state and appropriate statistical evidence linking the locus to sleep and frailty. A patient-level predictor needs a separate approved cohort with genetic and health outcomes. Training more epochs cannot replace those data.

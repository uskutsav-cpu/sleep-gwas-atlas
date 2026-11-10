# Small pretrained DNA model experiment

Written before fitting this new model. This experiment is retrospective in the broader project and is not a clinical sleep–frailty prediction study.

## Intended biological task

Predict whether a reference/alternate DNA pair was called statistically significant for an allelic effect in the Salomon human excitatory-neuron reporter experiment. Train on measured reporter calls, not AlphaGenome outputs, GWAS p values, source model predictions, gene labels, or annotation features derived from the measured outcomes. Nonsignificant calls are an imperfect comparison class, not proven biologically inactive DNA. The assay measures engineered reporter constructs, not endogenous regulation or patient health. No claim that stronger prediction is guaranteed.

## Sources and compact model

Source assay table: the archived `work/salomon_neuronal_variant_readout.tsv.gz`, originally from kircherlab/80K-Analysis at commit a686f1d552f97ebcc0224057b5a835002f7a3dc8. Only ref/alt DNA and measured assay labels are training inputs. Source software is MIT; underlying data terms and selection by prior sequence predictors are separately documented. Restrict this run to local research; do not redistribute raw assay sequences or imply a new unrestricted dataset license.

Base model: LongSafari/hyenadna-tiny-1k-seqlen-hf at commit e8c1effa8673814e257e627d2e1eda9ea5a373f6, BSD-3-Clause model card, roughly 0.45 million stored scalars. This public pretrained model learned genomic sequence patterns. The checkpoint includes human reference chromosomes later held out for the supervised assay task, so supervised chromosome holdout is not exclusion from unsupervised genome pretraining. Use a reviewed MLX port for local training with backpropagation; do not execute unreviewed remote model code. Verify numerical operations and gradients, document that official PyTorch executable parity is not established if that environment remains unavailable.

Implementation amendment before fitting: both attempted MLX runtime extractions failed for insufficient disk space and were removed. Use the much smaller MIT-licensed HIPS Autograd 1.8.0 package with NumPy automatic differentiation and its supported FFT derivatives instead. The same pretrained architecture and last-block fine-tuning task are retained; only the numerical backend changes. This is a NumPy/Autograd port, not an official HyenaDNA runtime. Record gradient checks and this implementation qualification.

Final backend amendment before fitting: replacing 74 byte-identical temporary NHANES source copies with hash-verified links to existing originals released sufficient space for MLX 0.22.1, which initialized the local GPU successfully. Keep the Autograd implementation as a numerical reference. The installed MLX real-FFT autodiff failed an identity-gradient check; use full complex FFT and real inverse output instead, which passed that check. The convolution and end-to-end parameter gradients must match the independent reference and finite differences before training.

## Data eligibility and split

- Keep autosomal, uniquely represented genomic contrasts with source-confirmed readout and an unambiguous significant/nonsignificant flag.
- Require exactly 270 bases for each allele, A/C/G/T only, exactly one substitution, matching source position, stated REF/ALT and orientation. Exclude conflicting or duplicated variant/orientation rows rather than treating them as independent people or observations.
- Training: chromosomes 1–16. Validation: chromosomes 17–18. Final test: chromosomes 19–22. Test chromosomes must not influence checkpoint or hyperparameter selection. Prior project models examined different subsets of this source; retain that history rather than claiming a newly collected independent study.
- Check canonical reverse-complement sequence-pair duplicates across splits. Remove training/validation examples that exactly match a later split and report exclusions. Chromosome separation does not guarantee no homologous sequence or other biological dependence.

## Fixed comparisons

1. Constant training-prevalence baseline.
2. Regularized logistic model on fixed 1–3-mer sequence-context and absolute allele-difference features, standardized using training only.
3. Frozen pretrained HyenaDNA encoder and a small trainable head.
4. The same head plus fine-tuning of the last pretrained block, using a smaller learning rate for pretrained weights.

The neural pair representation includes reference/alternate context and magnitude of allele embedding differences, with both DNA strands considered through reverse-complement augmentation and inference averaging if runtime permits. Inference should be symmetric to exchanging allele order for this binary magnitude task. Lock the architecture and learning rates before evaluating final-test metrics. At most eight fine-tuning epochs; validation average precision selects the checkpoint, with early stopping. Keep training curves and parameter-change checks to substantiate actual optimization.

Primary metric: average precision on the final supervised chromosome holdout, compared with its detected-call prevalence. Also report AUROC, Brier score for uncalibrated output, and precision/recall at a validation-selected threshold. Use uncertainty that acknowledges region dependence; chromosome-level resampling over four held-out chromosomes is a limited sensitivity, not precise generalization certainty. All compared models use the same eligible final-test pairs. Do not keep redesigning the model in response to final-test results.

Overall model recommendation is determined by the highest validation average precision among the simple k-mer baseline, frozen encoder/head and last-block fine-tuned model, with no final-test input into that choice. Save all models and test results, including when the more complex model fails to improve on a simpler one. Candidate predictions from any model are exploratory; neither a good validation score nor this recommendation authorizes disease-risk use.

## Project use

If reference/alternate contexts can be verified, the trained model may generate exploratory regulatory-effect rankings for the project's 33 variants. These predictions have no source assay labels at those exact candidates and cannot be treated as validated disease effects, causal genes, cell specificity across tissues, absolute risks, or genetic-score validation. A clinically useful joint sleep–frailty predictor still requires independent approved genotype–phenotype records.

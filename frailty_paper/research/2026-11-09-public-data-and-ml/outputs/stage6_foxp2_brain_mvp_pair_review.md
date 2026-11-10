# FOXP2 excitatory-neuron expression and MVP insomnia: expanded pairwise analysis

## Result

The expanded pairwise analysis remains **inconclusive**. SingleBrain FOXP2 expression in excitatory neurons produced one SuSiE credible set, but MVP insomnia produced no reported credible set on the same 857 variants. Consequently a multi-signal coloc posterior for these two traits could not be computed. This is not evidence that the expression and insomnia signals are distinct, nor evidence that they share a causal variant. The prior disease–disease finding and the separately conducted brain–HFRS analysis are unchanged.

## What changed

The previous three-way alignment required a variant to occur in both disease datasets as well as the expression data and the public reference. Here the comparison requires only the brain expression study, MVP and the NFE reference, while retaining the same 1-Mb GRCh37 region and non-palindromic SNP restriction. The pairwise input had 858 SNPs, 857 had exact position/allele matches in the public reference, and all 857 passed the pre-specified 0.20 allele-frequency-disagreement filter. These included all 720 SNPs of the preceding three-way analysis and rs12536335. The mean absolute MVP-versus-reference allele-frequency difference was 0.0050; this does not establish ancestry matching or eliminate LD differences. The matrix was complete, symmetric, unit-diagonal, and positive definite (minimum eigenvalue 3.11e-6). The custom gnomAD index decoder has not been independently verified with Hail.

The SingleBrain source lists insertion chr7:114534306:T:TA, but its fixed-effect beta and standard error both equal zero. It therefore has no usable fixed-effect estimate here and was not included; this is not a null expression result.

SuSiE-RSS used L=10, 95% credible-set coverage, fixed residual variance, the signed gnomAD NFE LD matrix, and effective MVP sample size 253,769. The expression sample-size sensitivity used 757 unique donors and 983 samples; the effective sample size is not established, so neither number should be interpreted as exact. Both fits converged for both sample-size settings. The brain fit returned one credible set (maximum PIP 0.602 and 0.614, respectively); MVP returned zero (maximum PIP 0.479). A maximum PIP does not by itself constitute an accepted credible set.

The strongest marginal expression variant among these 857 was rs1450833 (GRCh38 position 114225680; expression z=7.88, insomnia z=-1.78, both oriented to the MVP effect allele). The strongest marginal insomnia variant was rs2894699 (GRCh38 position 114419101; insomnia z=5.50, expression z=-4.17). These are descriptive marginal associations; differing lead variants alone do not prove distinct causal variants. The data are drawn from different samples and external LD; sample overlap, imputation, selection and the restricted variant set limit interpretation. Neither expression colocalization nor biological mechanism is established.

## Reproducibility

- Aligned data: `stage6_foxp2_brain_mvp_nfe_aligned.csv`; diagnostics: `stage6_foxp2_brain_mvp_pair_diagnostics.csv`; brain credible-set membership: `stage6_foxp2_brain_mvp_pair_members.csv`.
- Code: `stage6_foxp2_brain_mvp_pair_prepare.py`, `stage6_foxp2_brain_mvp_nfe_index.py`, `stage6_foxp2_brain_mvp_nfe_ld.py`, and `stage6_foxp2_brain_mvp_pair_coloc.R`. The index and LD receipt files record public source block hashes; the entire remote matrix was not downloaded or checksum-verified.
- Sources: [SingleBrain](https://singlebrain.org/), [gnomAD v2.1.1 LD resources](https://github.com/mkanai/ldcov), and the public MVP insomnia GWAS GCST90475826.

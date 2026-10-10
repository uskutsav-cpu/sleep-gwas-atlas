# Independent mathematical problem and identifiability review

2026-10-10. Scope: the new one-paper scientific steer, mathematical identification, and a bounded primary-source review. No simulations, fitted contrasts, estimators, source/reference/runtime census, individual data, or V4 edits were performed. The companion data contract and source ledger are part of this review.

**Decision: NO-GO for substantive method invention from the four proposed problems as presently specified.** Sampling-covariance propagation is established; unrestricted latent measurement error and transportability are unidentified; weak-denominator inference is a genuine limitation but no distinct unsolved contribution has been demonstrated. Conditional GO for designing a small, source-defined biological comparison using established inference, subject to the data and scientific gates below. This is neither an execution admission nor a Nature Communications submission recommendation.

**1. The strongest identifiable question.** For two oriented sleep measurements A and B and a fixed disease definition D, compare common-SNP genetic correlations

\[
\Delta=r_{AD}-r_{BD},\qquad
r_{iD}=c_{iD}/\sqrt{v_i v_D}.
\]

Here S contains additive common-SNP genetic variances v and covariances c for the declared phenotype transformations, covariates, source populations, and variant/LD target. With a common target genotype covariance R, an idealized definition is S_ij=β_i′Rβ_j. Population GWAS associations need not equal direct biological effects. Across different cohorts, interpreting all entries as one population S requires explicit compatibility/transport assumptions; source-defined effect similarities are the narrower target. Clinical insomnia, questionnaire trouble sleeping, usual duration, and rank-normalized actigraphy duration are different variables, not automatically exchangeable measures of one latent trait.

Partial participant overlap primarily changes estimation dependence. Under the cross-trait LDSC model,

\[
E[z_{ij}z_{Dj}]=\sqrt{N_iN_D}\,c_{iD}\ell_j/M
+N_{i\cap D}\rho_{iD}/\sqrt{N_iN_D}+a_{iD}.
\]

Separating slope from intercept requires the polygenic/LD model and adequate control of confounding; this is not a universal overlap correction. The original method treats overlap through the intercept and supports case-control rg without a liability-scale distinction under its model. [Bulik-Sullivan et al. (2015)](https://pmc.ncbi.nlm.nih.gov/articles/4797329/).

**2. Exact uncertainty problem, including h².** In a coherent three-trait joint fit, order θ=(v_A,v_B,v_D,c_AD,c_BD). For positive variances,

\[
g=\nabla\Delta=
\left(-\frac{r_{AD}}{2v_A},\ \frac{r_{BD}}{2v_B},\
-\frac{\Delta}{2v_D},\ \frac1{\sqrt{v_Av_D}},\
-\frac1{\sqrt{v_Bv_D}}\right)^\top.
\]

With V=Cov(θ̂), first-order uncertainty is

\[
\operatorname{Var}(\widehat\Delta)\simeq g^\top Vg
=\sum_i g_i^2V_{ii}+2\sum_{i<j}g_ig_jV_{ij}.
\]

Thus marginal rg SEs, two covariance SEs, or an overlap intercept do not supply all required cross terms. In correlation coordinates the same statement is Var(r̂_AD−r̂_BD)=Var(r̂_AD)+Var(r̂_BD)−2Cov(r̂_AD,r̂_BD). The shared D denominator's derivative vanishes at an equal-rg null, but the other h² derivatives and cross terms remain. A shared disease filename does not imply a single fitted v_D: original pair fits can use different masks and denominator estimates. Retaining those estimators instead requires a larger joint vector with both disease-denominator estimates and their covariance. Do not splice standalone h² estimates into pair-fit ratios or mix observed covariance with liability h². Freeze whether the reported statistic is the full-point ratio or a jackknife bias-corrected statistic.

Multivariable LDSC/GenomicSEM already estimates S and its joint sampling V; defined-parameter inference is an established delta/sandwich problem. The official tutorial explicitly warns that separate original LDSC outputs cannot populate off-diagonal V. [Grotzinger et al. (2019)](https://www.nature.com/articles/s41562-019-0566-x), [official genome-wide tutorial](https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Genome%E2%80%90wide-Models).

There is an important implementation distinction. Current `ldsc.R` forms `V_Stand` by congruent rescaling with point-estimated diagonal S. This holds denominators fixed rather than applying the complete ratio Jacobian. Its current code also uses retained-row quantiles separately for trait/pair blocks before stacking pseudovalues. Package output names therefore do not prove exact shared genomic deletions. [Official source, lines 243–244, 372–373, 430–441, 477–489](https://raw.githubusercontent.com/GenomicSEM/GenomicSEM/master/R/ldsc.R).

An explicit contrast defined from **unstandardized S,V** can propagate all denominator derivatives. Current DWLS `usermodel.R` computes sandwich uncertainty and a Jacobian for `:=` parameters. Its automated PSD/nonfinite fallbacks do not make an invalid V scientifically valid; they must not waive admission failures. The standardized-output branch is not the proposed full-denominator solution. [Official source, lines 156–217](https://raw.githubusercontent.com/GenomicSEM/GenomicSEM/master/R/usermodel.R). No package was installed or run. These mutable-source observations must be version-pinned before future execution.

**3. What an aligned jackknife would establish.** For common deletion events b=1,…,B, let θ_(−b) contain all components after the same genomic event. Under the usual equal-group/block-jackknife approximation,

\[
\widehat V=\frac{B-1}{B}\sum_b
(\theta_{(-b)}-\bar\theta_{(-)})
(\theta_{(-b)}-\bar\theta_{(-)})^\top.
\]

Alternatively recompute both nonlinear ratios inside each shared deletion and apply the formula directly to their difference. Finite-block direct nonlinear and first-order Jacobian results need not be identical. This is established resampling, not method novelty. The protocol must define weights, masks, intercept stages, deletion-specific versus fixed normalization, block balance, and the target of the jackknife. A common partition alone does not guarantee independent blocks or correctly estimated sampling uncertainty. Different masks may retain different SNPs inside a fixed block, but the deletion event's identity must persist. A common-SNP restriction changes the experiment and must be labeled.

Also distinguish participant-resampling uncertainty conditional on fixed effects/LD from uncertainty under a genomic random-effects model. By total variance, random-effect inference can include both expected conditional sampling variance and variation of the conditional mean across genomic architectures. Genomic blocks can reflect heterogeneous true signal as well as correlated estimation noise. The desired V and calibration experiment must name that repeated-sampling target; zero known participant overlap alone does not justify zero block-jackknife covariance.

Arrays numbered 1–200 do not identify those events. Permuting one vector preserves every marginal SE while changing its cross-covariance with another; the arrays alone cannot resolve the permutation. Reference-only partitions, invented controls, finite-panel decoder checks, and controller success supply no actual joint GWAS covariance. An empirical B-block V has rank at most B−1, favoring a prespecified small system rather than a 45- or 158-trait 200-block matrix. No optional V4 control is a prerequisite merely because it exists; choose the smallest scientifically justified established route once an eligible primary question is selected.

**4. Reliability does not identify a new rg correction.** Let Y_m=a_mT+e_m, a_m>0, with error genetically orthogonal to T, D, and the retained genotypes. Then

\[
\operatorname{Cov}_g(Y_m,D)=a_m\operatorname{Cov}_g(T,D),\qquad
\operatorname{Var}_g(Y_m)=a_m^2\operatorname{Var}_g(T),
\]

so r_g(Y_m,D)=r_g(T,D). After standardizing Y_m, the same positive scaling cancels. Classical error reduces observed-scale h² and precision, not population rg. Dividing rg by square-root phenotypic reliability is consequently incorrect under this model. Positive scalar effect-size dilution also cancels. Nonlinear transforms, thresholding, differential misclassification, heritable reporting behavior, and selection need not satisfy this model. Primary work already compares GWAS measurement quality and uses correlated jackknife uncertainty for repeated-versus-single-measure h². [Self-report inaccuracy study](https://www.nature.com/articles/s41562-024-02061-w).

The contrasting unrestricted model U_m=λ_mT+η_m permits genetically structured errors. It is not identified from S alone. Exact counterexample: let independent unit-variance genetic variables L,Q give D=L, A=L+Q, B=L−Q. The observed S is identical in all of these decompositions: T=L with errors Q,−Q (latent rg=1); T=L+Q with errors 0,−2Q (latent rg=1/√2); T=Q with errors L,L−2Q (latent rg=0). Every measurement loading is one. Nothing in these unchanged observed genetic variables tells us which component is biological sleep or measurement contamination. Equivalent unknown mixing persists even with noiseless SNP effect vectors.

Identification can be imposed by known calibration/anchors, restricted residual covariances, sufficiently informative repeated indicators, and fixed scale/orientation. Such assumptions are substantive; fitting a genetic factor does not validate them. Three nondegenerate indicators can identify a restricted one-factor genetic model, but an exactly identified fit cannot test its residual-independence assumptions. Two indicators need additional constraints or anchor/cross-covariance restrictions. Phenotypic test–retest reliability does not automatically identify genetically structured error, and repeated sleep observations also contain genuine temporal variation. PheMED already addresses relative effective dilution under effect-sharing assumptions; it does not identify unrestricted latent error or its biological cause. Its author documentation currently cites a preprint. [PheMED primary preprint](https://www.medrxiv.org/content/10.1101/2023.01.17.23284670v2), [author software](https://github.com/voloudakislab/PheMED).

**5. Ascertainment and transport are separate identification questions.** With a valid liability-threshold/ascertainment model,

\[
h_L^2=h_O^2\frac{K^2(1-K)^2}{p(1-p)\phi(t)^2},\qquad
t=\Phi^{-1}(1-K).
\]

K is population prevalence and p the sample-case fraction under the declared N convention. Coherent positive scale conversions give S*=DSD and leave rg unchanged. Unknown prevalence therefore need not prevent rg estimation under this model; it does prevent a justified absolute liability-h² interpretation. Effective-N inputs and sample prevalence must use a consistent convention. Congruent positive rescaling cannot repair indefinite S. The conversion does not remove unknown clinical selection or differential misclassification. [Lee et al. (2011)](https://pubmed.ncbi.nlm.nih.gov/21376301/).

If A is measured only in cohort 1 and B only in cohort 2, a contrast cannot separate measurement from cohort effects. For example observed correlations 0.6 and 0.2 fit either a −0.4 measurement effect with no cohort effect or a −0.4 cohort effect with no measurement effect. The missing crossed cells distinguish these worlds. A bridge with matched definitions across cohorts and/or multiple measurements within cohorts, followed by justified invariance assumptions, is needed. Matching ancestry labels alone does not provide it. Equality of selected-source rg is a transport description; population transport additionally needs a target population and identifiable selection/effect-modification model.

Similarly the overlap intercept identifies neither overlapping people nor their number: with N_1=N_2=100,000, overlap 20,000 with phenotypic correlation 0.5 and overlap 50,000 with correlation 0.2 both contribute 0.1. Shared confounding adds further ambiguity. Unknown overlap can be modeled for uncertainty without being biologically or administratively identified. Source independence remains a provenance requirement for replication.

**6. Weak h² is a real limitation, not an identified novelty opportunity.** For independent unit genetic Z_0,Z_1,Z_2, define D=Z_0 and U_i=√ε[r_iZ_0+√(1−r_i²)Z_i]. Then v_A=v_B=ε, c_iD=r_i√ε, c_AB=εr_Ar_B. As ε→0 every such S approaches diag(0,0,1), while Δ=r_A−r_B can differ throughout [−2,2]. At zero genetic variance rg is undefined. With nonvanishing estimation noise, arbitrarily sharp uniformly valid inference cannot recover the contrast near this limit. Clipping h², dropping negative deletions, selecting favorable high-signal realizations, or repairing eigenvalues does not solve identification/calibration.

Away from zero, delta inference is regular. Near the boundary, general profile/test inversion or projection of a valid joint confidence set onto Δ are established alternatives; sets may be broad or noninformative. Ordinary two-variable Fieller formulas are not automatically exact for c/√(v_iv_D), a nonlinear correlated denominator. A confidence-set projection inherits coverage only if its joint input set is valid; unknown V and finite LD blocks do not disappear. [Fieller/generalized ratio-confidence geometry](https://www3.stat.sinica.edu.tw/sstest/j19n3/j19n312/j19n312.html). SUPERGNOVA already explicitly discusses noisy heritability denominators; HDL/HDL-L supply likelihood approaches and LAVA supports local multivariable relations. None identifies latent biological measurement error. [SUPERGNOVA](https://link.springer.com/article/10.1186/s13059-021-02478-w), [HDL](https://www.nature.com/articles/s41588-020-0653-y), [HDL-L author documentation](https://github.com/YuyingLi-X/HDL-L), [LAVA primary paper](https://pubmed.ncbi.nlm.nih.gov/35288712/).

**7. Available-data and claim boundaries.** The dated V4 checkpoint records 190 reproduced native jobs, 158 standalone h² estimates, and 1,637 rg estimates. Their arithmetic is verified; it does not supply cross-fit covariance, latent measurement anchors, or individual overlap. The older 45-trait S/V and fixed-scale export are explicitly inadmissible for unqualified new SEM. Historical 41 validation estimates reuse sleep inputs; zero establish fully independent two-trait replication. FinnGen insomnia source-feasibility h² is not a new sleep–disease rg or a new source-admission decision. Active raw-to-munged work is preserved; its state is not rechecked here.

A supported new claim could be: a prespecified source/definition-specific rg contrast is inconsistent with equality, with calibrated joint uncertainty and independent validation of the appropriately narrow target. A nonsignificant contrast is not evidence of equivalence; equivalence requires a scientific margin and precision adequate to contain the interval within it. Neither conclusion identifies mechanisms, causes, latent reliability, or universal transport. Existing atlas associations are historical. Selecting pairs because either old P was favorable remains retrospective unless selection is appropriately addressed. Freeze a small finite family, source definitions, negative controls, precision requirements, and validation class before new outcomes; do not reuse old BH families as if they controlled a new contrast family.

**8. Adversarial decision.** The best argument for innovation is reliable inference for weak, heterogeneous denominator estimates with uncertain joint V. The counterargument is decisive at present: general ratio/test-inversion machinery already exists; alignment is an implementation/calibration requirement; and the harder biological targets need additional information rather than an estimator. A new contribution would first have to specify an identifiable sampling model, a demonstrated limitation of strong established baselines, and an improvement beyond standard delta, joint jackknife, or constrained confidence-set inference. This review provides no such gap. Do not start a new method or expensive simulation campaign now. Continue essential V4 completion and primary-source/phenotype/validation feasibility work; use the strongest established inference for a qualifying biological question. The one-paper goal remains possible in principle, but its new biological contribution and external validation are not established by this mathematical review.

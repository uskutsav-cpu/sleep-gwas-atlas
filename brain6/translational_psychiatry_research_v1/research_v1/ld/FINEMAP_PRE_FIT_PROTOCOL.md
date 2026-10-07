# New reference-qualified native fine-mapping protocol

Locked 2026-10-07 before any new GWAS/LD consistency diagnostic, SuSiE fit,
credible set or multi-signal colocalization outcome. Exact source archives
are being reacquired under the parent's identity/permission manifests. This
is a retrospective exploratory method-validation experiment using the same
discovery studies, not independent GWAS replication or prospective discovery.

## Scientific target and primary units

Question: within the previously reported chr5 insomnia–ADHD interval, does a
source-qualified multiple-signal model support one or several shared effects,
and do alternative priors/model sizes materially alter that interpretation?
Window fixed GRCh37 chr5:103447968–104447968. This is candidate C; its rs2431108
lead strongly tags published rs77960 in new reference LD (already observed),
so no first-discovery or distinct-signal conclusion follows from lead naming.

Primary units: two trait fits on the identical full-window intersection and
all coloc-SuSiE signal-pair comparisons. No smaller P-selected window, SNP
significance filter, outcome-selected signal count or favorable comparison
selection. No confirmatory family/P threshold applies to these exploratory
posteriors; the separate four-hypothesis independent-validation FWER family
remains 0.05/4. New chr6/chr11 reference resources do not admit disease fits
under this protocol without a distinct source-qualified locked experiment.

## Exact model and source gates

Pinned native susieR0.14.2 and coloc5.2.3, R4.6.1; source tarballs and every
required dependency recorded in the CRAN manifests. Use original discovery
insomnia and Demontis2023 ADHD archives with exact historical SHA256 identity.
Respect ADHD source restriction against public result-file redistribution;
raw regional inputs and any reconstituted source association table stay local.
New analysis aggregates/posteriors are derived research output.

Require original source-supported chr/position/allele conventions, beta and
SE on the named effect allele, effect/model semantics, and per-variant N.
For ADHD beta=log(OR), SE=reported logOR SE, n_eff=4/(1/Nca+1/Nco), control
A1 frequency for QC. Insomnia model and reported N must be verified from
original source README/publication; if case/control effective N cannot be
constructed from the supplied literal analyzed N and documented case/control
ratio, no N-assisted RSS fit is admitted. An n=Inf z-only sensitivity cannot
launder that source uncertainty into a primary result.

Allele-harmonize to GRCh37 REF/ALT 1000Gv5b503EUR reference keys. Exact rsID
identity additionally checked for named leads from primary Ensembl mappings;
v5b removed other rsIDs, so full identity uses coordinate AND allele pair,
not coordinate alone. Reject palindromes, non-biallelic SNPs, duplicates,
contradictory keys, invalid/nonfinite beta/SE/P/N, INFO<0.9 where reported,
reference/source MAF<0.01 where reported, and absolute source-vs-reference
ALT-frequency difference>0.15 where frequency exists. Never infer missing
source frequency or INFO; report that trait's unverified QC field explicitly.
No outlier/SNP can be removed merely to increase a posterior.

Require >=500 shared admitted SNPs and >=80% of eligible nonpalindromic,
MAF-qualified reference SNPs represented by each GWAS; >=80% of each trait's
valid nonpalindromic biallelic SNPs in the window must overlap the reference.
The union/intersection denominators are retained. Require >=80% of otherwise
eligible SNPs within +/-20% of median n_eff; excluded N outliers and all
reasons are preserved. This is a single-N RSS approximation: even passage
does not establish an exact common-sample regression likelihood. Require
both traits min P<5e-8 in the full admitted window to interpret shared
credible signals. Absence fails adequacy, not a substantive biological null.

Use actual float64 signed ALT dosage Pearson LD, no clipping, shrinkage or
PSD projection. Require reference protocol QC passage. Rank deficiency with
p>503 is expected and permitted by the PSD RSS method; do not pretend this
changes the earlier frozen Cholesky-required 0/32 failure.

## Pre-fit consistency diagnostic and stop criteria

Run susieR::estimate_s_rss(null-mle) using median effective N plus n=Inf
sensitivity, recording exact LD/GWAS hashes. Run kriging_rss with that s and
retain every SNP diagnostic. These are model diagnostics, not GWAS results.
Require s<=0.10 and no suspected allele-switch marker satisfying documented
logLR>2 and |z|>2. If flagged, stop this native posterior branch and report
LD/summary-statistic incompatibility. Only an externally documented source
or harmonization error permits correction; require a new versioned input
receipt and regression test. Do not correct a SNP from diagnostic likelihood
alone. Record whether finite reference sampling or true allele mismatch is
unresolved. This conservative operational hold is not a universal SuSiE
calibration guarantee. Reference503 may inadequately capture full source LD;
scalar/matrix numerical agreement cannot establish summary-statistic fidelity.

## Fits, convergence, priors and fixed sensitivities

Primary susie_rss(z=ALT-aligned beta/SE,R=actualLD,n=median n_eff,L=10,
estimate_residual_variance=FALSE,coverage=.95,min_abs_corr=.5,
max_iter=1000,tol=1e-3,refine=FALSE). Full covariance untouched. Save full
native R objects locally, all variant PIPs and all credible sets, purity,
ELBO history, convergence and warnings. Nonconvergence, ELBO drops>1e-6,
nonfinite probabilities or absent qualifying credible sets are explicit
failures/unresolved, not successful null fine-mapping.

coloc.susie primary p1=p2=1e-4,p12=1e-5. Report every signal pair, H0–H4,
including H3 vsH4. Fixed p12 sensitivities1e-6 and5e-5. Descriptive support
requires each trait convergence and >=1 qualifying95%CS,purity>=.5,
H4>=.8 and H4/(H3+H4)>=.8; no posterior declares a causal gene/variant.
Run fixed L=1 andL=5 fits plus n=Inf sensitivity for both traits, same
variants/reference, max_iter/tol/convergence gates, all signal-pair outputs.
Run coloc.abf on original beta/SE across the same shared SNP universe with
same three priors as an explicit single-causal-variant sensitivity. These
comparisons assess model dependence; they are neither independent cohorts
nor independently calibrated significance tests.

Document exact CS memberships and whether named prior lead markers appear;
high LD cannot discriminate their causality. Discovery cohort overlap,
binary-GWAS approximation, source/reference sample-size difference,
post-selection and missing independent two-trait data remain limitations.
No historical output, canonical LAVA status or ABF label is changed.

# Executed reference-qualified chr5 fine-mapping

Candidate C now has native multiple-signal-capable analysis, using restored
original insomnia/ADHD sources and newly reconstructed signed genotype LD.
This materially improves method interpretability over the inherited
single-signal ABF analysis. It strengthens a previously published shared
region; it does not establish a newly discovered locus, independent
replication, causal variant, effector gene or psychiatric mechanism.

## Completed empirical computations

- Three true-genotype regional reference reconstructions: 503 EUR samples
  and 2694/4878/4627 admitted variants on chr5/6/11, respectively.
- Two raw source-qualified chr5 input/harmonization audits, retaining 2185
  common variants after prespecified INFO/frequency/N/identity rules.
- Two native RSS consistency diagnostics plus both large-N sensitivities.
- Eight native SuSiE-RSS fits: primary L10, fixed L1/L5 and L10 large-N,
  for each of insomnia and ADHD. All converged without ELBO decreases.
- Twelve native coloc-SuSiE invocations (four fixed model settings by three
  priors), each producing one retained signal-pair comparison.
- Three native single-signal ABF sensitivities on that identical SNP universe.
- Twelve independent Python recomputations of H0–H4 from native signal Bayes
  factors, maximum absolute posterior disagreement 3.442e-15.
- Six EUR/subpopulation rs2431108–rs77960 reference comparisons and 24
  additional Lin-prior/reference comparisons; these are LD diagnostics.

Software/API regression tests are counted separately: ten source-free tests
and three native tests using synthetic data. They do not validate a biological
finding. Three full LD matrices independently rebuilt from compact genotype
factors match their original SHA256 exactly.

## Primary result

|Quantity|Insomnia|ADHD|
|---|---:|---:|
|Common variants|2185|2185|
|Median effective N supplied|312492.707847208|128213.795046423|
|RSS consistency s|0.0140860905313902|0.00917219053539645|
|Potential allele-switch flags|0|0|
|Converged primary iterations|2|6|
|Qualifying 95% credible sets|1|1|
|Primary credible-set size|70|8|
|Minimum absolute CS LD correlation|0.646449155326643|0.727916896727039|
|Maximum variant PIP|0.137300679996284|0.2947218560012|

The ADHD credible set is entirely contained in the insomnia credible set.
Both include rs2431108, published rs77960 and published insomnia–MDD rs30266.
The exact membership table remains authoritative; rs171697 is palindromic
and was excluded from primary trait harmonization. Its reference-only dosage
LD was evaluated separately after a labeled retrospective protocol extension.

Primary coloc-SuSiE H4 is 0.991568229456347 at p12=1e-5; H3 is
0.00842968755924881. At the lower p12=1e-6, H4 is 0.921629512291253 and H3
0.0783511270652316. At p12=5e-5, H4 is 0.998302193479262. Every fixed L1,
L5 and large-N sensitivity also exceeds the descriptive H4 and H4/(H3+H4)
thresholds. Across model settings, lower-prior H4 ranges from
0.919491100946998 to 0.921675808038437. Same-universe ABF H4 is
0.991868276548707 by default and 0.924228165890470 at the lower prior.

No second qualifying credible signal was detected with L up to 10 in this
admitted window. That is a model-conditional finding; it cannot rule out
poorly tagged, rare, low-powered, excluded palindromic or external-window
causal effects. Low maximum PIPs and broad shared LD leave individual causal
variants unresolved. Both primary fits are approximate binary-GWAS RSS.

## Source and model qualifications

Exact archive identities were recovered: insomnia SHA256
32848cd92a6324c9048cad4a804cf1546299af1d0281981c6ca43c2046a9065b;
ADHD SHA256 c58a96031ec44b1edba81c91d603037a2efd03fcd946531424866e3f5f1d40c5.
These are the same discovery sources, not new psychiatric cohorts. Original
raw row counts differ from historical filtered counts; they are not treated
as source identity failures when exact bytes/SHA match.

Insomnia beta is log(OR) for source A1. Its literal per-variant analyzed N
was transformed by the documented UKB case fraction 109402/386533 to an
approximate effective N. Per-variant case fraction is not supplied. ADHD
uses source A1 logOR/SE and variant Nca/Nco-derived effective N; contributing
meta-analysis cohorts and participant counts vary. Reference ancestry is
European, but 503 reference individuals do not recreate UKB or the ADHD
meta-analysis population exactly. Reference sampling, residual cohort
structure, unknown cross-study overlap and case/control model approximation
still limit posterior calibration.

Source OR and SE rounding leads to source-P versus beta/SE-Z discrepancies
(maximum absolute Z discrepancy 0.09582 insomnia and 0.03740 ADHD). The
prespecified analysis retains literal beta=log(OR)/SE-derived signed Z and
source P for association admission, rather than silently replacing effects
with P-derived statistics. Numeric compatibility is adequate under the
operational gates; the gates do not prove exact model correctness.

Inherited canonical LAVA failure and original 0/32 UKB-factor LD failures
remain unchanged. The new independent genotype reference and distinct native
experiment do not retrospectively repair those historical gates. Chr6/chr11
reference availability alone does not admit disease fine-mapping; insomnia
source association minima in those A/B windows do not meet the component
association requirement, and their source/model questions remain separate.

## Independent automated adversarial conclusions

Different lead SNP names cannot support chr5 novelty: rs2431108 and Zu's
rs77960 have EUR r²=0.9954806027688181. Lin's rs171697 has r²=
0.9955156955525478 with rs2431108; Lin's rs30266 is perfectly correlated
with rs77960 in EUR. These are reference relationships, not conditional
causal-signal tests. Native H4 support is prior-stable within the fixed
sensitivity grid but remains conditional on the GWAS/LD model and marker
universe. A separate read-only automated review recomputed all8native-fit
PIP/alpha probability and ELBO checks, all8credible-set memberships/minimal
95%coverage and purity, and all12posterior vectors using an explicit
i≠j Bayes-factor grid. The shared2185-SNP matrix exactly matches the full
2694-SNP submatrix. Maximum purity disagreement was6.78e-15; posterior
agreement was3.45e-15. Its authoritative record is
`../statistics/INDEPENDENT_NATIVE_FINE_MAPPING_REVIEW.md`. This is a
separate automated implementation/review, not independent human validation
or empirical posterior calibration. One detected credible signal per trait
is not independent replication. No human statistical-genetics validation
was performed.

Before promotion, an independent human reviewer should evaluate source N
semantics, case/control RSS approximation, finite-reference LD uncertainty,
source OR/SE quantization, missing case-fraction variation, unknown overlap
and full-window completeness. Truly independent two-trait replication,
a separate larger eligible reference and independently supported molecular
colocalization would change the scientific interpretation most decisively.

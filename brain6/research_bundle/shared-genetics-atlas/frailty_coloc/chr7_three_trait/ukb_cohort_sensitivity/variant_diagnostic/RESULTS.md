# What drives the FI-source sensitivity at chromosome 7?

A post hoc decomposition of the [UKB-only comparison](../RESULTS.md) finds that the two FI inputs give similar leading candidate sets, while the FI–ADHD shared-signal posterior depends on small tail weights outside those sets. These are conditional, single-signal approximate-Bayes-factor weights within the 1,055 selected variants; they are not independently validated causal probabilities.

| Diagnostic | Meta-analysis FI | UKB-only FI |
|---|---:|---:|
| FI 95% conditional candidate-set size | 25 | 26 |
| FI candidates also in ADHD's 34-variant 95% set | 0 | 0 |
| Normalized FI–ADHD variant-weight overlap | 0.0000462 | 0.0005843 |
| Conditional H4/H3 odds at p12=10⁻⁵ | 0.0462 | 0.5847 |
| Leading variant conditional on FI–ADHD H4 | rs9969232 | rs9969232 |
| Its conditional share of the H4 weight | 0.216 | 0.260 |

All 25 meta-analysis FI candidates occur in the UKB-only FI set of 26. FI candidate sets overlap the insomnia set (25 and 26 variants respectively); the analogous normalized FI–insomnia weight overlap changes only from 0.02111 to 0.02055. These numerical intersections depend on the 95% cutoff and do not demonstrate that the effects are identical, distinct, or restricted to those variants. In particular, zero intersection does **not** imply zero overlap of the full weight distributions or rule out a shared signal.

For two traits, the normalized weight overlap is the sum of products of their conditional single-signal weights. With this model's p1=p2=10⁻⁴ and p12=10⁻⁵, multiplying overlap/(1−overlap) by 1,000 reproduces the H4/H3 odds. Thus a roughly 12.7-fold change in a very small weight overlap helps explain the posterior shift; it does not determine why the input weights differ. The overlap is not an LD correlation, a probability of biological interaction, or evidence for a shared causal gene. All calculations hold the selected support fixed; alternative boundaries and omitted or multiple signals remain untested.

The [calculation](analyze.py) and [machine-readable output](diagnostic.json) reproduce the parent posterior odds and verify normalized weights. The matched input remains on the SSD under the ADHD source's terms. Source sample composition, statistical precision, study dependence, allele conventions and FI sleep items remain important limitations; this decomposition does not resolve them or failed LAVA QC.

An [input-swap diagnostic](counterfactual/RESULTS.md) shows how this particular posterior reacts to FI beta and standard-error columns. Its hybrid inputs are artificial and have no sampling interpretation; they cannot identify why the source fits differ.

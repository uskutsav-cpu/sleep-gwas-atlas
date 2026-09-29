# Exploratory GWAS ABF fine-mapping and trait–trait coloc

The frozen primary LAVA family did not pass. None of these 25 pair-specific
PLACO candidates has been promoted. These results are a separate,
prospectively locked **single-causal-variant ABF** analysis and cannot be
described as LAVA confirmation, multi-signal SuSiE fine-mapping, or
independent replication.

## Input and method decision

The initial SuSiE-RSS plan in `brain6_exploratory_finemap_coloc_v1` admitted
27 trait-locus inputs and five shared-pair inputs, but **0/32** regional LAVA
UKB reference matrices passed the pre-specified numerical LD gate. Raw
absolute correlations reached 1.01956. Representative raw minimum
eigenvalues were −0.05069 (chr2, 2,092 SNPs), −0.01660 (chr11, 1,493 SNPs),
and −0.14594 (chr17, 3,592 SNPs). No SuSiE PIP, SuSiE credible set, or
SuSiE-coloc posterior was calculated. A 503-person 1000 Genomes EUR PLINK
panel independently available in the project has only 219–585 eligible
common non-palindromic SNPs per protected region (21/25 below 500 before
GWAS matching), so it cannot supply dense multi-signal LD for this branch.

The separate v2 ABF method and priors were committed before association
posterior inspection. ABF uses source beta and SE and does not require a
regional LD matrix, but it assumes at most one causal variant per trait in
each interval. Its 50 trait-locus units produced 28 ABF evaluations:
24 model-conditional 95% credible sets and four with association posterior
below the frozen 0.80 credible-set gate. Nine units failed dense input QC and
13 Dashti long-sleep units remain blocked by model/effect-scale semantics.
The 24 credible sets range from 3 to 1,445 variants (median 53); a large set
does not fine-map a variant precisely. The full 56,826 variant-level ABF PIPs
are checksum-bound in the SSD detail file named by `integrity_provenance.json`.

## Trait–trait coloc

Six insomnia–ADHD candidate intervals passed both-trait dense input QC and
have full H0–H4 coloc.abf posteriors. All other 19 pair units have explicit
NOT_RUN reasons: six lack dense MDD input and 13 involve the unresolved
long-sleep model. The table below reports the frozen p12=1e-5 analysis and
the prespecified p12=1e-6 sensitivity; the complete H0–H4 vectors and
p12=5e-5 sensitivity are in `trait_coloc_25.tsv` and
`trait_coloc_prior_sensitivity.tsv`.

| Protected interval (GRCh37) | Shared SNPs | PP.H3 | PP.H4 | PP.H4 at p12=1e-6 | Frozen descriptive ABF label |
|---|---:|---:|---:|---:|---|
| chr11:28.08–29.08 Mb | 1,579 | 0.7186 | 0.0029 | 0.0003 | Not met |
| chr14:98.14–99.14 Mb | 3,108 | 0.2956 | 0.5454 | 0.1071 | Not met |
| chr3:117.14–118.14 Mb | 2,340 | 0.0779 | 0.9207 | 0.5373 | ABF model support; prior-sensitive |
| chr5:103.45–104.45 Mb | 2,359 | 0.0089 | 0.9911 | 0.9174 | ABF model support across tested priors |
| chr5:87.19–88.19 Mb | 1,420 | 0.5946 | 0.1379 | 0.0157 | Not met |
| chr7:113.57–114.57 Mb | 1,155 | 0.1472 | 0.8432 | 0.3496 | ABF model support; prior-sensitive |

Three of six meet the **predeclared ABF-model descriptive** threshold at
p12=1e-5, but only the chr5:103.45–104.45 Mb result retains PP.H4>=0.80
under the low p12 sensitivity. None has undergone valid multi-signal
conditioning or formally confirmed local shared-genetic validation. The
posteriors depend on the assumed one-causal-variant model and prior; no gene,
causal variant, or disease mechanism is established from these values.

## Files and audit

`fine_mapping_50.tsv` has exactly one row per candidate and trait;
`fine_mapping_credible_sets.tsv` contains every reported 95% ABF credible-set
member. `trait_coloc_25.tsv` has one row per protected pair candidate,
including PP.H0–H4 only for the six executed pairs. The independent
`validate_brain6_exploratory_abf_v2.py` audit checks all 50/25 identities,
source input hashes, allele matching, posterior normalization, conditional
PIP and credible-set rules, all prior-sensitivity vectors, and reproduces the
ABF null and H0–H4 calculations independently of the R package. Its passing
checksum record is `integrity_provenance.json`.

The primary method source is [Giambartolomei et al. (2014)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383).
The limitation from multiple causal variants is documented in
[Wallace (2021)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009440).

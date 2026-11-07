# Four-trait marginal threshold screen

This post-hoc marker screen compares FinnGen R13 registry insomnia, FI, FFS and FinnGen R12 hospital frailty on the same 1,138,262-marker sleep/FI/FFS HapMap3 support. The pre-result [protocol](PROTOCOL.md) retains 4,553,048 marker–trait slots in the Bonferroni denominator, giving p=1.098×10⁻⁸. HFRS had unique allele-compatible values at 1,137,927 markers; 335 had no HFRS match and remain unknown. All 504 reported threshold-crossing markers have unique allele-compatible FI positions, HFRS matches and matching chromosome numbers. Coordinate/allele joining does not replace a full reference-aware cross-build audit.

| Observed marginal threshold pattern | SNPs |
|---|---:|
| FFS alone | 288 |
| HFRS alone | 171 |
| Insomnia alone | 23 |
| FI alone | 19 |
| FI and FFS | 3 |
| Insomnia and any frailty phenotype | 0 |

These are **SNPs, not independent loci**, and “alone” describes a significance pattern rather than proven trait specificity. Sample size, phenotype, ascertainment, alleles and LD affect which p-values cross. HFRS and FinnGen insomnia have overlapping participants; FI and FFS share UK Biobank. The first three p-values are conservative bounds based on the rounded aligned Z scores; HFRS uses displayed source p-values with finite precision. No marker meeting an individual threshold for both insomnia and frailty does not exclude shared genetic architecture, nor does a marginal pattern prove a shared causal variant. This screen cannot identify HFRS-only variants outside the three-trait common support. It is separate from PLACO+ and enrichment-based conjFDR.

The smallest HFRS p-value among the 171 HFRS-only threshold patterns is 3.64×10⁻³⁰ for rs4420638, GRCh37 chr19:45422946. Its other p-values are approximately 0.00285 for insomnia, 0.988 for FI, and 0.0441 for FFS; none passes the four-trait cutoff. This is an example, not a newly established HFRS locus or evidence of no effect on other traits. Neighboring and correlated threshold-crossing markers require appropriate LD analysis and independent source validation before locus-level claims.

`crossings.tsv` retains effect alleles, Z scores, source HFRS effect and p-value, aligned effect, and FI-derived position for every crossing; `summary.json` records hashes and coverage. The three-trait component of the selected set exactly matches the earlier three-trait table re-filtered at this stricter threshold (333 markers); `verification.json` records that comparison and checks the HFRS-only rows. The earlier [selected 368-marker HFRS lookup](../marginal_priority/hfrs_followup/RESULTS.md) did not contain the additional HFRS-only crossings. Neither analysis evaluates participant-level item-excluded phenotypes.

An independent threshold-first join checked the HFRS source separately and reproduced all 171 HFRS crossings on the common support. Another ten HFRS source SNP IDs meet the same numerical cutoff but are absent from the three-trait support; they are listed in `independent_hfrs_check.json`. They are not part of this four-trait family and the lookup does not establish their joint trait patterns or a valid correction for a broader HFRS scan.

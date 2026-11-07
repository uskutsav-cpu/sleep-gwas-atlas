# Post-hoc sensitivity to extreme LD–Z flags

Removing the **10 distinct variants** with |standardized kriging residual| >5 in any of the three original eQTL, insomnia or ADHD diagnostics left **1,538 of 1,548** variants. This common exclusion removed neither rs62264764 nor the dominant adrenal eQTL marker rs6804394. The threshold and this comparison were chosen after reviewing the pilot diagnostics. The flagged variants are not established errors, and their removal is not a correction or a formal calibrated test.

| Result | Original support | After removal |
|---|---:|---:|
| Adrenal LSAMP credible-set size | 3 | 3 |
| Insomnia credible-set size | 75 | 75 |
| eQTL–insomnia SuSiE H4 at p12=10⁻⁵ | 0.001430 | 0.001431 |
| eQTL–insomnia SuSiE H4 at p12=10⁻⁶ | 0.00014314 | 0.00014331 |
| ADHD qualifying credible sets, median N | 0 | 0 |
| ADHD qualifying credible sets, minimum N | 0 | 0 |

The two retained credible sets had exactly the same members as before. Their SuSiE fits converged, and the same-support single-signal results were similarly close; the full [fit summary](fit_summary.tsv), [credible sets](credible_sets.tsv), [coloc status](coloc_status.tsv), [posterior table](coloc_insomnia_1e-05.tsv), and [same-support ABF results](single_signal.tsv) are available. The comparison indicates limited sensitivity to *this particular removal rule*, not robustness to all LD mismatches, missing variants, sample-size assumptions, priors, or alternative biological signals. There is still no ADHD–eQTL coloc posterior from the pilot because the ADHD fit returned no qualifying set.

Recomputing the diagnostics after removal still produced one insomnia and two ADHD variants with |standardized residual| >5; their values and the remaining sub-threshold distribution can change when the conditioning set changes. Threshold counts are descriptive, not validated error rates or proof of allele flips. The strongest residual fell from 11.27 to 5.46 in ADHD, from 6.60 to 5.06 in insomnia, and from 5.40 to 4.34 in the eQTL. Reference LD remains out-of-sample and singular; the input study frequencies, missing actual per-variant model semantics, GTEx sample size, small reference sample and selected window retain all limitations described in the [original pilot](../RESULTS.md). The result remains **NOT ADMITTED** as validated colocalization and does not repair the LAVA QC gates.

[Protocol](PROTOCOL.md), [removed-variant diagnostics](removed_diagnostics.tsv), [retained variants](retained_variants.tsv), [rerun code](run.R), [diagnostic summary](kriging_summary.tsv), and [independent subset verification](verification.json) preserve the selection rule and data hashes. Posterior normalization and set membership were checked; they do not validate the model. No AlphaGenome outputs were used to fit these models.

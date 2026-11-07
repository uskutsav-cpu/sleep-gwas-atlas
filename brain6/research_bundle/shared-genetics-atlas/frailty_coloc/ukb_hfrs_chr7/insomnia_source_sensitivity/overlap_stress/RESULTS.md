# Correlated-statistic stress: UKB HFRS and questionnaire insomnia

The selected chromosome 7 two-trait HyPrColoc result is **sensitive to the assumed correlation between association statistics**. On the same 1,125 variants, the model returned an HFRS–insomnia cluster at an assumed effective correlation of 0, but no two-trait cluster at the tested positive values of 0.05 and above. The true correlation is unknown, so this does not establish which result is correct or that the traits do not share a signal.

| Hypothetical effective correlation | Returned cluster | Model posterior | Candidate SNP |
| ---: | --- | ---: | --- |
| −0.3 | HFRS, insomnia | 0.9831 | rs2894699 |
| −0.1 | HFRS, insomnia | 0.9550 | rs2894699 |
| −0.001 | HFRS, insomnia | 0.7159 | rs2894699 |
| 0, independence path | HFRS, insomnia | 0.7106 | rs2894699 |
| +0.001 | HFRS, insomnia | 0.7051 | rs2894699 |
| +0.05, +0.1, +0.2, +0.3, +0.5 | No two-trait cluster returned | Not reported | Not reported |

“No cluster” is the package's thresholded output under that particular hypothetical setting, not a negative biological result. The near-zero values approach the independent path, a useful software continuity check. Negative settings are mathematical stress cases, not estimates of the actual sign. We did not estimate a switch point or a correlation from the data.

The package multiplies `trait.cor` by `sample.overlap`; we set the latter to a matrix of ones and vary their working product. This does not assert complete overlap. The real number of shared participants and the relevant residual correlation are unverified, and the genome-wide LDSC intercept is not a validated substitute for this local parameter. The model uses a default identity LD matrix, one signal per trait, an unverified HFRS raw-effect scale, and fixed illustrative priors and thresholds. Its posterior and thresholded cluster output use a different conditional prior and method from the [earlier coloc H4](../RESULTS.md); the values should not be equated. Posterior values are not multiple-testing-corrected p-values.

Both input traits use UK Biobank, but exact cohort overlap is unknown. This test strengthens the caution around the earlier independence-based shared-signal hypothesis; it does not estimate a valid overlap correction, establish a causal variant or FOXP2 mechanism, rule out secondary signals, or repair the failed LAVA gate. See the [pre-run protocol](PROTOCOL.md). All ten fits returned with no recorded R warnings; the matched-input hash was verified, and an independent check confirmed near-zero continuity. These are software checks, not calibration.

The [all-scenario table](../../../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-overlap-stress-clusters-20261028.tsv`), [complete saved model outputs](../../../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-overlap-stress-20261028.rds`), [run log](../../../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-overlap-stress-20261028.log`), [script](../../../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-overlap-stress-20261028.R`) and [hash receipt](../../../../PUBLICATION.md#external-dependencies) (external target: `/Volumes/Extreme SSD/brain6-hfrs-insomnia-overlap-stress-receipt-20261028.json`) are on the SSD.

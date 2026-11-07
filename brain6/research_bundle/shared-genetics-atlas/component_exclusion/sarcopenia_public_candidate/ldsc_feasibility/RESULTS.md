# Sarcopenia-factor LDSC feasibility: not a validated item-excluded frailty result

Four LDSC fits returned numerical estimates on 984,470 common markers. They use the public sarcopenia-factor summary statistics and the existing questionnaire-insomnia, registry-insomnia, walking-pace and grip sources. **These are feasibility diagnostics, not calibrated estimates for physical frailty with sleep or exhaustion removed.** The factor's exact input GWAS, their coding, loadings, sample overlap and the meaning of its constant N=522,632 remain undocumented in the checked release.

| Comparison with released factor | Model-reported r_g (SE) | LDSC-reported p | Cross-trait intercept (SE) |
| --- | ---: | ---: | ---: |
| Questionnaire insomnia | −0.1688 (0.0207) | 3.24×10⁻¹⁶ | −0.0639 (0.0075) |
| Registry insomnia | −0.1289 (0.0287) | 7.26×10⁻⁶ | 0.0073 (0.0075) |
| Walking pace | +0.3369 (0.0190) | 3.31×10⁻⁷⁰ | 0.2753 (0.0101) |
| Left grip strength | +0.7749 (0.0084) | Printed as 0 | 0.5519 (0.0155) |

LDSC printed zero for the grip p-value due to numerical underflow; it is not literally zero. Multiplying the three nonzero p-values by four gives approximately 1.30×10⁻¹⁵, 2.91×10⁻⁵ and 1.33×10⁻⁶⁹. This numerical four-test benchmark is conditional on the unverified model and does not correct the much broader search, phenotype selection or prior analyses. It should **not** be used to declare discoveries. The two sleep estimates were not formally compared.

The positive walking and grip estimates provide numerical orientation context if the other source traits increase with better performance. They do not establish what higher values of the latent factor mean: the factor loadings and its input list are unavailable, and those physical traits may themselves overlap the factor's inputs. The large walking and grip intercepts emphasize the need to investigate dependence and model fit; an intercept is not a participant-overlap count. The questionnaire sleep and physical traits include UK Biobank; registry insomnia is Finnish, so phenotype differences, overlap, sampling, and the general-European LD reference all matter. The registry insomnia univariate intercept was 1.0765 (SE 0.0095).

The assumed constant factor N gave an observed-scale h² of 0.1817 (SE 0.0084) and intercept 0.9635 (0.0301) in each fit. Neither number validates the sample-size meaning, genetic factor or LDSC applicability. The [Genomic SEM methods documentation](https://github.com/GenomicSEM/GenomicSEM/wiki/4.-Common-Factor-GWAS) discusses how factor scaling and effective sample size affect heritability calculations; this release does not document its factor-specific scaling or N derivation. Do not use these estimates as proof of a sleep-independent frailty mechanism, causality, a locus, independent replication or resolution of failed LAVA QC.

## Checks and reproduction

Each of four complete source streams matched the deposited MD5, checked 282 concatenated gzip members, and retained 984,470 exact-coordinate, non-palindromic common markers with zero eligible duplicate IDs. Another 910 otherwise matching rows failed numeric/frequency filtering. All four archived input hashes matched their earlier receipts. The full source was streamed through an anonymous pipe; it was not saved. Source hashes and counts, each complete LDSC log, [protocol](PROTOCOL.md), [streaming code](stream_source.py), [runner](run.sh), [calculation and logs](summary.json), and [input verification](input_verification.json) are retained in this folder. The very first attempt failed before reading the source because of a local path error; the corrected `run2_` outputs are the only results used.

These checks establish data identity, matching and numerical completion. They do not validate the factor, its inputs or the statistical calibration. The requested same-sample sleep- and exhaustion-excluded FI/FFS GWAS remains outstanding. The [parent source review](../RESULTS.md) records release quality and selected-marker results.

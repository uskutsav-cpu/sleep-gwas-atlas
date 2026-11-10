# Reference-group sensitivity of the TPMI 17-marker screen

The conditional 17-marker weighted-Z screen was repeated with subsets of the same 1000 Genomes phase 3 EAS reference panel, keeping the TPMI browser summary results fixed. Each reference subset supplies both allele-dosage SDs and within-chromosome dosage correlations. Cross-chromosome correlations were set to zero, as in the original screen. All 17 markers varied in each subset.

| Reference samples | People | Insomnia 327.4 weighted Z | Organic/persistent insomnia 327.41 weighted Z |
|---|---:|---:|---:|
| All EAS | 504 | 2.44 | 1.28 |
| CHB + CHS | 208 | 2.48 | 1.27 |
| JPT | 104 | 2.51 | 1.27 |
| KHV + CDX | 192 | 2.36 | 1.28 |

CHB and CHS are the panel's Beijing and Southern Han Chinese samples; JPT is Japanese in Tokyo; KHV is Kinh in Ho Chi Minh City; CDX is Chinese Dai in Xishuangbanna. The subsets are components of the all-EAS panel, so the rows are not independent replications. The small reference samples are not TPMI participants or an exact target-specific LD estimate. Their similar estimates only indicate limited sensitivity to these particular reference substitutions; they do not establish statistical calibration, target generalizability or prediction. This calculation holds the target association results fixed and does not capture uncertainty from the target GWAS, allele coding, imputation, selection or sample overlap.

The [original screen](stage6_tpmi_17_marker_exploratory_z_review.md) remains exploratory and assumes browser beta is an ALT effect; five markers of the complete count remain unverified. The tests and outcomes are correlated and postselected. The [machine-readable output](stage6_tpmi_reference_population_sensitivity.json) and `python3 outputs/stage6_tpmi_reference_population_sensitivity.py` record the checks and reference groups; the script first rebuilds the original 17-marker result. TPMI requests acknowledgement and citation of [Chen et al.](https://doi.org/10.1038/s41586-025-09350-y) and [Yang et al.](https://doi.org/10.1038/s41586-025-09680-x).

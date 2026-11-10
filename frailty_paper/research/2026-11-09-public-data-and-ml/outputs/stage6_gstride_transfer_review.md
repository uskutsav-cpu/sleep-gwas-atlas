# Exploratory gait-feature transfer from FRAILPOL to GSTRIDE

[GSTRIDE](https://zenodo.org/records/6883292) is an openly licensed (CC BY 4.0) dataset of 163 older adults with foot-sensor gait features and Fried frailty items. It is distinct from the Polish [FRAILPOL](https://doi.org/10.6084/m9.figshare.c.7874411) dataset, although individual or relative overlap was not independently verified. GSTRIDE has no verified linked genome-wide genotypes or insomnia outcomes, so this analysis cannot validate the project's genetic scores.

To investigate the earlier model's generalizability, we fitted a **new**, simpler three-feature multinomial model on all 668 FRAILPOL participants and evaluated it on GSTRIDE. The inputs were average stride time, speed and stride length. FRAILPOL's left and right values were averaged; GSTRIDE's single-foot averages were used directly. Training-only means and standard deviations were applied to GSTRIDE. Its 0–5 Fried count was grouped as 0 robust, 1–2 prefrail and 3–5 frail. Feature definitions and cohorts were inspected before this exploratory test; it is not an untouched, preregistered evaluation of the earlier 26-feature model.

| Measure on 163 GSTRIDE participants | Result |
|---|---:|
| Accuracy | 0.387 |
| Balanced accuracy | 0.322 |
| Always predict the largest GSTRIDE class | 0.472 |
| Recall: robust / prefrail / frail | 0.954 / 0.013 / 0.000 |

The model predicted almost everyone as robust, and identified none of the 21 frail participants. This is a warning about this particular transfer, not a conclusive test of all gait models or proof of causation. The 0.472 majority benchmark uses the observed target class distribution and is descriptive, not a deployable rule learned in advance.

The features are **only approximately comparable**: FRAILPOL derives two-foot ankle features from a timed up-and-go task, while GSTRIDE uses one foot and its own walk, sensors and processing pipeline. Gait speed is also related to the frailty phenotype's slowness item. The countries, ages, settings, test protocols and class frequencies differ; diagnostic coding and sample independence cannot be assumed identical. Feature mismatch, domain shift, imbalance and post hoc selection can contribute to poor performance. We did not tune on the GSTRIDE data, but they have now been inspected and must not be reused as a pristine test. A separate verified independent cohort and standardized feature protocol would be required for an unbiased future assessment.

The 1.37-GB Zenodo ZIP was not downloaded in full. An HTTP byte-range ZIP reader extracted only the public `Database_register.csv` and gait README, with Python's ZIP CRC check on extraction. The repository advertises the full archive checksum, which was **not** verified because the archive was not downloaded. `Rscript work/frailpol_gstride_transfer.R` reproduces the model from locally extracted data; [metrics](stage6_gstride_transfer_metrics.csv) and [confusion counts](stage6_gstride_transfer_confusion.csv) are aggregate only. The original extraction scripts are `work/gstride_remotezip_inventory.py` and `work/gstride_extract_manifest.py`. No genetic or sleep data, individual predictions, or clinical advice were produced.

## Feature-distribution comparison

A descriptive follow-up compared the three modeled features before scaling. The GSTRIDE means were higher than the FRAILPOL means by **0.96 training standard deviations** for stride time, **0.64** for speed, and **1.23** for stride length. For example, average stride length was 0.616 m in FRAILPOL and 0.912 m in GSTRIDE. These are differences in the recorded features from two different tasks and processing pipelines; the measures have not been established as interchangeable. The comparisons describe overall samples, not matched individuals or group-adjusted differences, and do not identify whether measurement, population, frailty, or other factors caused the shifts or the model's poor performance. No correction or model tuning was performed based on them.

`Rscript work/frailpol_gstride_feature_shift.R` writes the [aggregate feature comparison](stage6_gstride_feature_shift.csv). It uses no participant predictions and does not change the earlier model.

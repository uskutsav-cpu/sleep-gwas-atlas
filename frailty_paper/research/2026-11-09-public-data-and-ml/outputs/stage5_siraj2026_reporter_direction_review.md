# Reference/alternate direction cross-check for five reporter candidates

For the five variants meeting the authors' selected-library active-and-skew rules in SK-N-SH, `stage5_siraj2026_five_variant_direction_check.csv` compares source `Log2Skew` with a separately computed aggregate ratio:

`log2((mean_RNA_alt / mean_Plasmid_alt) / (mean_RNA_ref / mean_Plasmid_ref))`.

All inputs are positive and come from the authors' selected-library GTEx processed table. Their [public analysis code](https://github.com/julirsch/finemapped_mpra/blob/main/code/preprocess/GTEx_munge.R) constructs the labeled mean reference and alternative RNA and plasmid columns from the oligo metadata and count data. For all five selected variant/library pairs, the sign of this aggregate ratio is the same as the reported `Log2Skew`, supporting an ALT-versus-REF orientation for these **five particular records**. The simple ratio is not a reproduction of the generalized linear model or its uncertainty and is not proof of the convention for every source row. All five records had previously passed reciprocal GRCh37/38 and reference-base checks; the source measurements use the GRCh37 labels.

| GRCh38 variant | Aggregate reporter ALT/REF direction | Predicted SK-N-SH DNase ALT–REF direction |
| --- | --- | --- |
| chr3:117873501:C>A | positive | negative |
| chr5:88401026:G>A | negative | negative |
| chr15:74179953:A>G | negative | negative |
| chr15:74193421:G>T | negative | negative |
| chr15:74214138:C>T | negative | negative |

The observed ratio and model DNase sign agree for four of these five records and differ for one. **This is a descriptive, postselected comparison**, not an estimate of prediction accuracy or independent validation: reporter RNA activity from an isolated construct and DNase accessibility are different measurements, and the model track and reporter measurements are not necessarily the same samples or conditions. The five were chosen precisely because their reporter allelic activity passed the study's criteria; they are not a representative benchmark or independent risk loci. Full model training and benchmark overlap has not been excluded, and ALL_FOLDS is not a held-out sequence model. Matching a sign does not identify a target gene or show a shared insomnia–frailty effect.

The SMAD5-region alternative chr5:136207208:G>A is outside the authors' selected analysis subset and is not among these five. No inference about its reporter allele direction follows from this check. Other cell lines, primary-cell measurements, and source records have not been checked by this count-based method.

Reproduction: the source processed table, checksum and selected-library extraction are documented in `stage5_siraj2026_best_library_source_receipt.json` and `stage5_siraj_all_cell_selected_extract.R`. Run `Rscript outputs/stage5_siraj_best_library_extract.R GTEX_RDS selected_with_counts.csv` and then `python3 outputs/stage5_siraj_reporter_allele_orientation.py --selected-counts selected_with_counts.csv`; the latter also requires `stage5_siraj2026_allcell_reporter_model_comparison.csv`. The underlying model track values are in the live AlphaGenome score archive.

A later check across all 284 selected cell-line result groups found nine cases where the aggregate count-ratio sign differs from the fitted skew, including rs60271 in HepG2. All nine were outside the study's active-and-skew calls, but the source of the discrepancy is unresolved. This reinforces that the count-ratio check cannot establish the GLM's allele convention generally. See `stage5_siraj2026_allcell_orientation_review.md`.

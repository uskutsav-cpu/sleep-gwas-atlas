# Public reporter validation across cell lines

We applied the published [Siraj et al. analysis code](https://github.com/julirsch/finemapped_mpra/blob/main/code/preprocess/read_in_functions.R) to the study's two processed data files: `best_library == 1`, excluding `other_test`, followed by the authors' active-element and allelic-skew rules. Exact GRCh37 REF/ALT matching was restricted to 433 queue variants with reciprocal, reference-checked mappings. The remaining 55 were not evaluated. These are externally generated reporter measurements, but independence from all AlphaGenome training and benchmark resources has not been established.

The two source subsets yielded 565 rows for **42 distinct queued variants**. Rows repeat across gene, tissue, trait and study annotations and are not independent assays. After collapsing identical result groups within source, cell type, library and reported values, `stage5_siraj2026_allcell_reporter_model_comparison.csv` contains 284 groups. These retain the number and identities of contributing source records, exact allele, reporter results, cell-line label, selected AlphaGenome DNase score, and project GWAS context. All groups had at least one model DNase track with the same visible cell-line label, although the experiment and model track are not the same samples or assays.

| Reporter cell line | Variants in selected subset | Variants meeting active-and-skew rules |
| --- | ---: | ---: |
| SK-N-SH neuroblastoma | 42 | 5 |
| K562 blood-derived line | 42 | 3 |
| HepG2 liver-derived line | 42 | 3 |
| HCT116 colon-derived line | 28 | 2 |
| A549 lung-derived line | 23 | 0 |
| GM12878 lymphoblastoid line | 23 | 0 |

There are five distinct variants meeting the rules in at least one line; the same variant can appear in multiple lines and analysis sources. Classification follows the provided source columns but retains the authors' best-library and source-context filters. A variant not classified as an emVar is not proven to have no effect; coverage, detection power, activity and cell context matter. These immortalized cell models do not replace normal brain or primary immune tissue.

The previously identified SMAD5-region alternative has a significant raw SK-N-SH library result but is **not in this selected analysis subset**; it should not be counted among the five. rs323509 has a selected SK-N-SH row but does not meet the reported activity-and-skew rules. The [neuronal review](stage5_siraj2026_neural_reporter_review.md) explains both in context.

`stage5_siraj2026_allcell_coverage.csv` records all 488 queue alleles by cell line, including mapping exclusions and lack of selected records. RNA reporter skew and AlphaGenome DNase are distinct molecular measurements; the original A/B allele convention is not verified globally. A separate count-based cross-check supports ALT-versus-REF orientation for the five selected SK-N-SH candidates only; see `stage5_siraj2026_reporter_direction_review.md`. Neither a matched cell-line name nor genomic proximity proves a target gene, tissue-specific mechanism, disease colocalization or independent model validation. Some source data and genomic sequences may overlap AlphaGenome's training or evaluation; `ALL_FOLDS` is not a held-out sequence model.

To reproduce, run `Rscript outputs/stage5_siraj_all_cell_selected_extract.R GTEX_RDS GTEX_CSV` and again with the trait RDS, using the source files and hashes in `stage5_siraj2026_best_library_source_receipt.json`. Then run `python3 outputs/stage5_siraj_allcell_integrate.py --gtex GTEX_CSV --traits TRAITS_CSV --store LIVE_SCORE_DIRECTORY`. `stage5_siraj2026_neural_reporter_review.md` describes retrieval of the raw source and mapping assumptions. The source rules were used descriptively after selecting the project's variants; no pre-registered or calibrated model-performance test is claimed.

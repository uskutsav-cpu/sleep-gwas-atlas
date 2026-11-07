# Comparing AlphaGenome and genetic rankings

The selected [reporter benchmark readiness audit](reporter_benchmark_readiness/RESULTS.md) finds no positive rows under a conservative activity-plus-skew proxy and no validated negatives. A classification performance claim is therefore not supported by that subset.

The 425 variants were selected from six Brain6 single-signal credible-set unions before the predictions were reviewed. This is a descriptive within-region comparison, not a genome-wide benchmark or a test of model performance. It compares three different quantities:

- **PIP:** a trait-specific, conditional single-signal statistical fine-mapping probability under the fitted model; it depends on the GWAS and region model.
- **AVI:** a broad overall AlphaGenome variant-impact score from a third-party mirror, not a tissue-specific disease score. Values and allele-order encoding have been spot-checked, not independently verified for every variant.
- **Frontal-cortex DNase change:** the absolute value of log₂[(ALT+1)/(REF+1)] for the summed predicted signal within 501 bp. It is a local, tissue-track prediction and not an experimentally measured effect, p-value or direct gene-level score. The signed values compare genomic ALT with REF, not risk direction.

| Candidate | Region's selected variants | Absolute DNase rank; signed ratio | AVI raw-score rank | Insomnia PIP rank | ADHD PIP rank |
|---|---:|---:|---:|---:|---:|
| rs62264764, chr3 | 128 | 1; −0.2734 | 1 | 6 | 2 |
| rs2431108, chr5:103.45–104.45 Mb | 71 | 69; +0.00053 | 1 | 1 | 3 |
| rs1476535, chr7 | 66 | 11; −0.0207 | 14 | 1 | 9 |
| rs4916723, chr5:87.19–88.19 Mb | 3 | 2; −0.0124 | 1 | not estimated | 1 |

The rs2431108 example illustrates why the scores should not be treated as interchangeable: it ranks highly in overall AVI and in both genetic PIP lists but has a small predicted change in the chosen adult frontal-cortex DNase window. This does **not** imply that AVI is wrong or that the variant has no effect in another molecular track or tissue. rs62264764 combines the strongest local DNase prediction in its region with high ADHD PIP ranking, but is sixth by insomnia PIP and its region's shared-signal colocalization is prior-sensitive. Neither pattern establishes causality.

Of the top-five ranked SNPs per score, chr3's DNase list shares one SNP with the ADHD PIP list (rs62264764) and none with the insomnia PIP list; the chr5:103 Mb DNase list shares none with either PIP top five. On chr5:87 Mb, all three available SNPs are necessarily in every top-five list, so the apparent overlap is uninformative. Other region/trait results and boundary ties are in `overlap_summary.tsv`. These data were selected using genetic evidence, and ranks can reflect LD, region choice and ties; no p-values, random-background expectations, or multiple-testing corrections are claimed. Missing PIPs are excluded from their rank denominators rather than set to zero, and trait rankings are withheld where the original association-posterior gate did not admit a credible set.

`variant_ranks.tsv` gives each variant's values and ranks, `compare.py` reproduces the tables, and `input_hashes.json` binds the selected prediction and region files. This analysis does not train an ML model. The user's approved non-commercial AlphaGenome access terms apply to the outputs.

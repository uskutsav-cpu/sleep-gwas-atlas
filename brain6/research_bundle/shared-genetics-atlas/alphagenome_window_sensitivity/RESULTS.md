# AlphaGenome accessibility summary: window and offset sensitivity

The prior pilot summarized each saved predicted adult frontal-cortex DNase track as log₂[(ALT sum + 1)/(REF sum + 1)] over 501 bases centered on the variant. This descriptive sensitivity analysis recalculates that same statistic for 101, 501 and 2,001 bases and offsets of 0.1, 1 and 10, using the existing predictions for all 425 selected variants. These choices were made after seeing the original results; they are not independent predictions or a prespecified robustness test. All comparisons are among this selected set, not genome-wide, and the offsets and windows change the meaning and scale of the statistic.

| Centered window | Offset | rs62264764 score | Rank by absolute score among 425 |
|---:|---:|---:|---:|
| 101 bases | 0.1 / 1 / 10 | −0.498 / −0.387 / −0.120 | 2 / 1 / 1 |
| 501 bases | 0.1 / 1 / 10 | −0.292 / −0.273 / −0.167 | 1 / 1 / 1 |
| 2,001 bases | 0.1 / 1 / 10 | −0.109 / −0.107 / −0.0918 | 1 / 1 / 1 |

At the 101-base, 0.1-offset setting, rs4372430 has the largest absolute score (+0.513). rs62264764 remains negative and ranks first or second across these nine settings. At offset 1, however, 66 of the 425 variants change sign between the 101- and 501-base summaries, and 78 change sign between the 2,001- and 501-base summaries. These include potentially tiny numerical differences: sign changes are not tests of meaningful biological effects. In particular, choosing a window or offset after inspecting the outcomes can favor a desired ranking.

The outputs are model-predicted aggregate track signal differences, not calibrated variant-effect probabilities, measured chromatin accessibility, evidence of a gene target, or proof of a shared disease mechanism. A negative ALT-versus-REF summary does not indicate a protective or risk direction for insomnia or ADHD. Larger windows can dilute or mix local positive and negative changes. The original 1,048,576-base model inputs are unchanged; only saved 2,001-base output windows were resummarized. The JBrowse service's model version is not exposed, limiting reproducibility across future service changes.

`run.py` verifies each saved compressed track hash and confirms it reproduces the original 501-base, offset-1 statistic; `window_scores.tsv` retains the raw REF/ALT sums and all nine summaries. This is descriptive analysis of predictions already obtained under the user's non-commercial authorization, not training another model. The [pilot report](../../brain6-alphagenome-pilot/results/jbrowse_expanded/EXPANDED_REPORT.md) documents selection and track provenance.

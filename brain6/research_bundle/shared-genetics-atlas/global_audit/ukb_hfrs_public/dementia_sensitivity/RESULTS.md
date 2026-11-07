# UK Biobank HFRS with dementia weights removed

Using the authors' [public summary statistics](https://doi.org/10.6084/m9.figshare.31259023.v1), we compared the original hospital frailty risk score (HFRS) with the authors' version excluding dementia weights. Both use UK Biobank data and the same **1,020,621 markers**. This is a dementia-component sensitivity, not the requested FI without sleep or Fried score without exhaustion.

| Insomnia source | Original HFRS r_g (SE) | HFRS without dementia r_g (SE) | Original minus modified (paired SE) | Contrast p, Bonferroni across two |
| --- | ---: | ---: | ---: | ---: |
| Questionnaire | 0.5288 (0.0334) | 0.5281 (0.0323) | 0.0007 (0.00334) | 1.000 |
| Registry | 0.3755 (0.0416) | 0.3718 (0.0410) | 0.0037 (0.00336) | 0.543 |

Both modified-score correlations remain positive and meet Bonferroni correction across the four individual tests. Neither paired difference meets the separate two-contrast correction. The corresponding **unadjusted** 95% intervals for original minus modified are −0.00585 to 0.00725 for questionnaire insomnia and −0.00289 to 0.01029 for registry insomnia. These results suggest little change in these particular genome-wide point estimates after removing the dementia weights; they do not prove zero effect, rule out changes at specific loci, or show a sleep-independent frailty mechanism. The differences use the rounded LDSC log estimates; their uncertainty uses a paired 200-block jackknife.

This is same-data sensitivity, not independent replication. The questionnaire insomnia and HFRS sources share UK Biobank, with exact overlap unverified. Registry insomnia comes from FinnGen, and a general-European LD reference may not adequately reflect Finnish ancestry. Correlation is not causation, a shared causal variant or a local-correlation result. The adjusted p-values apply only to the stated post-hoc families, not the wider atlas.

## Source and reproducibility

The complete 280,547,193-byte modified-source gzip was streamed and matched its publisher MD5 `68d413550cbc2b3d7b4857ad36b7b4f7`; SHA-256 `8174fac2f512cf8129f1f745fe4d12f087f0e792083d1b5b20a15037a570c9c5`. All 1,020,621 selected original-HFRS markers matched by ID, reference position and allele pair and passed the specified quality filters. No matching IDs were duplicated; all retained BETA/SE-derived normal p-values agreed within 0.02 on the log10 scale with the source. These checks do not validate the upstream model or prove allele semantics. All four fits retained the same markers, shared-insomnia heritability delete vectors matched exactly within each paired contrast, and individual jackknife SEs matched the LDSC logs to displayed precision.

All 9,812,854 rows in the modified source, like the original file, report **N=407,643**, whereas the authors' article and release description state **407,463**. The reason remains unknown; the analysis uses the file's value, and source-level N is not independently verified. A first attempt to save the subset ran out of space and left an incomplete file on the small SSD workspace; a separate valid filtered file was saved at `/Volumes/Extreme SSD/brain6-ukb-hfrs-without-dementia-20261028.sumstats.gz`. No partial output was used or deleted.

The [protocol](PROTOCOL.md), [source-audit code](stream.py), [source receipt](stream_receipt.json), [run script](run.sh), [paired calculation](summarize.py), [numerical results](results.json), and per-pair logs and delete vectors preserve the analysis. The original HFRS inputs and logs are documented in the [parent analysis](../RESULTS.md). Matched sleep- and exhaustion-excluded FI/FFS GWAS and a valid local-correlation family remain outstanding.

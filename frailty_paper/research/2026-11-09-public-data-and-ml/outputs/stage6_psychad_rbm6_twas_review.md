# Exploratory RBM6 brain mural-cell association

## Findings

A European-ancestry PsychAD prediction model for genetically regulated RBM6 expression in brain mural cells met the source study's expression-model filter (31 SNPs; cross-validated R²=0.0362; prediction q=2.71e-7). Using the model's own covariance matrix and the project's existing public GWAS statistics, S-PrediXcan gave the following **exploratory, uncalibrated** results. Negative Z means lower model-predicted expression corresponds to higher GWAS trait values under the aligned alleles; it does not demonstrate that changing RBM6 expression causes the trait.

| GWAS | Model SNPs used | Z | Two-sided p |
|---|---:|---:|---:|
| Atkins frailty index | 28/31 | -3.761 | 0.000170 |
| Mak frailty index | 31/31 | -5.585 | 2.33e-8 |
| Mak modified physical frailty | 31/31 | -4.226 | 2.38e-5 |
| Ye Fried frailty score | 31/31 | -3.769 | 0.000164 |
| Jansen insomnia symptoms | 28/31 | -2.247 | 0.0246 |
| FinnGen R13 clinical insomnia | 30/31 | -0.455 | 0.649 |
| FinnGen R12 hospital frailty score | 30/31 | -1.693 | 0.0905 |

The stronger FI and physical-frailty associations are **gene-model associations**, not individual-level validation. The insomnia estimates differ in strength, and neither FinnGen result reaches nominal p<0.05; those results are not proof of no effect. This model and these genes were selected after examining project evidence. The seven tests, earlier analyses, correlated traits, sample overlap and model selection mean these raw p values cannot be treated as a prespecified discovery or as independent replications. UK Biobank contributes to multiple inputs; FinnGen samples may overlap between its insomnia and hospital-frailty studies. The 3p21 region contains multiple genes and correlated variants, so a TWAS signal may reflect local LD or other mechanisms rather than mediation by RBM6. Expression-model training overlap with other data was not established. The hospital score and physical/deficit frailty definitions differ; sleep-related items or codes can also overlap the outcomes.

## Sensitivity to model variants

Using the 28 SNPs shared by all seven inputs preserved the direction of all estimates. The most influential leave-one-out removal was rs2353586: the Mak FI Z changed from -5.585 to -3.203 and the Jansen insomnia Z from -2.247 to -0.939. Removing a SNP changes the prediction model; this is not a causal or conditional analysis. Details and all results: `stage6_psychad_rbm6_sensitivity_review.md`.

## Inputs and quality checks

The model and covariance are from the [public PsychAD snTIM release](https://zenodo.org/records/17435917); source model quality criteria and limitations are in the [study](https://www.nature.com/articles/s41586-026-10836-6). ZIP-member CRCs and hashes were checked on targeted downloads, but not the complete multi-gigabyte archive. The RBM6 covariance has 496 upper-triangular entries for the 31 model variants. The reference model represents approximately 779 donors and is specific to European-ancestry brain mural cells; this is not a test of all cell types or ancestries.

GWAS matching required an exact rsID, unique GRCh37-to-GRCh38 positive-strand mapping when relevant, model coordinate and allele-pair agreement, a non-palindromic SNP, finite beta and positive SE. Inputs were oriented to the model's effect allele. The association cache was extracted earlier using the project's regional windows; matching and exclusions appear in `stage6_psychad_rbm6_gwas_harmonization.csv`. The full FinnGen hospital-frailty file was re-scanned and its published MD5 verified; its receipt is `stage6_psychad_hfrs_extract_receipt.json`. GWAS effect-allele definitions and provenance inherit the earlier project source audits.

Official [MetaXcan](https://github.com/hakyimlab/MetaXcan) at commit `e069063a10539fe92adadb645fc667a40c2cf885` was run on a one-gene copy of the original model. The model's `gene_name` and `n_snps_in_model` fields were copied into names expected by MetaXcan. One compatibility change makes a one-element NumPy array explicitly scalar, recorded in `stage6_metaxcan_numpy_compatibility.patch`; it does not alter the intended formula. No GWAS heritability or sample-size calibration was supplied: MetaXcan explicitly warns its resulting Z and p values are uncalibrated for inflation. Independent NumPy calculation of the published summary-association formula matched the program within 6.2e-8 Z units; the 31-SNP covariance minimum eigenvalue was -8.3e-17, consistent with numerical roundoff. This arithmetic check does not validate source assumptions, population matching, sample overlap or causal interpretation.

Machine-readable results: `stage6_psychad_rbm6_twas_exploratory.csv`; formula check: `stage6_psychad_rbm6_formula_check.csv`; program logs: `stage6_psychad_rbm6_metaxcan_run_receipt.json`. Reproduction scripts are `stage6_psychad_extract_small_models.py`, `stage6_psychad_extract_mural_cov.py`, `stage6_psychad_rbm6_prepare_inputs.py`, `stage6_psychad_hfrs_extract.py`, `stage6_psychad_rbm6_harmonize.py`, `stage6_psychad_rbm6_run_metaxcan.py`, and `stage6_psychad_rbm6_formula_check.py`. The harmonizer requires `work/stage5_env/bin/python` with pyliftover and the previously cached region data. Further testing requires a frozen gene/cell/trait set, appropriate correction and calibration, independent outcome cohorts, and robust multi-signal colocalization before any shared-risk or mechanistic claim.

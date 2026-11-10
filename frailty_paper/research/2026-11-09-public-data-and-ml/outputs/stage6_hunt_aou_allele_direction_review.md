# HUNT–All of Us insomnia candidate direction: source-code audit

29 October 2026. Source-code review supports, but does not fully authenticate every deployed source table, the assumption that the All by All EUR single-variant beta is for the API's ALT allele. Under that assumption, the published HUNT insomnia OR and All of Us EUR insomnia beta have the same sign on the HUNT effect allele for five of six selected candidates. The two with unadjusted insomnia p < 0.05 in both cohorts—rs113851554 and rs4886140—are among those five. These are postselected variant associations; agreement alone is not independent score validation, proof of shared causality or a formal replication result.

## Trace of the reported direction

- The All by All authors' [v8 SAIGE loader](https://github.com/atgu/aou_gwas/blob/8b71fb31bdb3d570334cef8ebf807c4e4d4838ba/v8/pipeline/saige_aou_v8.py#L582-L590) stores its Hail allele key as `[Allele1, Allele2]` and carries the `BETA` field into the variant table. [SAIGE's own documentation](https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE) states that its single-variant BETA refers to Allele2.
- The [deployed browser release's slow interval handler](https://github.com/broadinstitute/all-by-all-aou-browser/blob/317fe0c9369aeef30c6a64c35d9f64e00c3c613a/axaou-server/src/variants/annotations.rs#L620-L654) reads the ancestry/phenotype Hail table and returns its typed beta and allele fields. Its [pinned genohype decoder](https://github.com/broadinstitute/genohype/blob/3b901ea564fdb463e307f95232b6479cc166c349/core/src/genomic/extract.rs#L102-L135) maps the first and second allele-array elements to REF and ALT, and reads the `BETA` value. The browser code commit matches the public release identifier displayed by the browser. This is a source-code trace; it does not independently inspect the protected deployed Hail table schema or exclude intervening source changes.
- The independent [assembly and pair check](stage6_lane_aou_six_liftover_review.md) maps all six assumed GRCh37 Lane positions to the reported GRCh38 positions on the forward strand and finds matching base pairs. It does not resolve every multiallelic representation or sample relationship.

| Variant | HUNT effect allele | HUNT OR (95% CI) | All of Us ALT | All of Us beta as returned | Sign on HUNT allele, assuming ALT coding |
|---|---|---|---|---:|---|
| rs113851554 | T | 1.13 (1.06–1.20) | T | +0.124754 | same |
| rs4886140 | G | 1.04 (1.01–1.07) | G | +0.0421471 | same |
| rs1942262 | A | 1.05 (1.01–1.08) | A | +0.0218737 | same |
| rs9845387 | C | 1.09 (1.02–1.16) | A | −0.0315585 | same after allele swap |
| rs10865954 | T | 1.01 (0.98–1.04) | C | −0.0181298 | same after allele swap |
| rs35881094 | G | 1.02 (0.99–1.05) | G | −0.00491398 | opposite |

The HUNT ORs are rounded; several intervals include 1 and several All of Us p-values are nonsignificant. A five-of-six sign count is not an independent inferential test, and reversing an effect allele requires the exact biallelic contrast to be verified. rs113851554 and other project loci have multiallelic representations; do not use this audit as a blanket authorization for score-dosage flipping. Study populations, ascertainment and relatedness may differ or overlap. The [browser walkthrough](https://allbyall.researchallofus.org/walkthrough) cautions that beta magnitude is not readily interpretable given the association-testing weights; do not equate its displayed value with HUNT's OR or use it as a validated score weight. Variant associations do not establish a frailty association or clinical prediction.

The [calculation table](stage6_hunt_aou_allele_direction_sensitivity.csv) includes p-values and the explicit signed transformation. Run `python3 work/cross_cohort/hunt_aou_direction_check.py` to reproduce it from existing extracts; [checksums](stage6_hunt_aou_allele_direction_receipt.json) identify the inputs and output. Source code was inspected at the pinned revisions linked above.

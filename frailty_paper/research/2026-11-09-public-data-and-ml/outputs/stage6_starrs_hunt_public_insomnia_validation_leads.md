# Public insomnia evidence search: STARRS and HUNT

Checked 29 October 2026. Neither source below is an independent validation of the project's 15- or 22-marker score.

## STARRS

[Stein et al.](https://doi.org/10.1038/s41380-018-0033-5) studied lifetime insomnia disorder in three US Army STARRS samples (3,237 cases and 14,414 controls overall) across European, African and Latino ancestry groups. This is a potentially useful *different* sleep phenotype, but the young military sample and case-enrichment strategy limit generalization to frail older adults. Participant/relative overlap with MVP or other sources requires checking; do not assume it is zero.

We searched all rsIDs in the article's six main-text table elements in its [Europe PMC full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6129221/fullTextXML) against the exact rsIDs in both current prototype files. None matched. These are selected published table entries, not a genome-wide inventory; zero matches do not show lack of association, coverage, or absence from supplementary results. A complete, accessible STARRS insomnia summary-statistics file was not located in this search. The [Partners study](https://doi.org/10.1038/s41598-020-63792-0) reports using STARRS summary results, but this is not an independent download link or evidence of permission for this project. The STARRS study reports an association for a **different** genome-wide UK Biobank insomnia PRS; it is not evidence for either project score, FI or FFS.

## HUNT

[Lane et al.](https://doi.org/10.1038/s41588-019-0361-7) report testing a different 57-variant sleep score in HUNT (14,923 cases and 47,610 controls) and Partners. This establishes a prior comparator rather than validating our score. HUNT has a [public summary-statistics repository](https://dataverse.no/) for some studies, but a keyword search of its API on 29 October did not establish a complete HUNT-only insomnia GWAS file. Search indexing and terminology can miss material, so absence from this search is not proof that no file exists. If an appropriate full file is located, verify provenance, access conditions, allele/build, age and phenotype, sample overlap and complete candidate coverage before any fixed-score-related association lookup. Summary associations cannot assess person-level discrimination, calibration or incremental prediction.

## Reproducibility

Input prototype files: `unvalidated_15_variant_score_prototype.csv` and `unvalidated_22_variant_score_prototype.csv`. Source XML saved as `work/starrs_insomnia.xml`; the HUNT keyword-query response saved as `work/hunt_insomnia_dataverse_search.json`. No human data were downloaded and no contact was made.

## HUNT4 physical-frailty and FI opportunity, with fee restriction

A [HUNT4 70+ frailty study](https://doi.org/10.1007/s40520-024-02839-y) reports 9,956 participants assessed using both Fried criteria and a 35-item FI. HUNT has also supported insomnia-phenotype and genome-wide genetic research, including Lane's score study cited above. These separate publications make HUNT scientifically attractive for comparing FI and physical frailty, but their sample sizes are not interchangeable: the intersection of HUNT4 70+ frailty, relevant sleep measures and imputed genotypes is unverified. FI sleep-item overlap, variant coverage, discovery overlap and target-population eligibility also need checks.

The [HUNT Data Bank 2026 price list](https://www.ntnu.no/documents/10304/1300106453/Prisliste%2BHUNT%2Bdatabank%2Bog%2Bforvaltning_2025.pdf/f3ac8eb0-7bc2-28b2-32c4-d01be6cfdc71?t=1743582519139) specifies fees that differ by funding/project category, including charges for provision of large GWAS datasets and small SNP/PGS datasets. It is not evidence that our institution is eligible for an exception or that costs would be zero; project-specific costs and conditions remain to be determined. Thus HUNT is not currently a verified free individual-level validation route and should not displace genuinely no-fee or open options. A public summary-statistics file for this frailty comparison has not been established.

The [HUNT4 item audit](stage6_hunt4_frailty_sleep_item_audit.md) confirms that the 35-item FI includes both insomnia and tiredness, and identifies a published correction for an earlier grip-strength data error. Use corrected values and item-excluded sensitivity analyses if this cohort ever becomes accessible.

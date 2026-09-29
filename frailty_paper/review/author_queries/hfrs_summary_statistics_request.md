# Access query draft: published FinnGen HFRS summary statistics

**Status:** Draft for human review; not sent. Intended recipient: FinnGen results-support contact or the corresponding author(s) of the HFRS GWAS.

**Study:** Mak et al., “Large-scale genome-wide analyses with proteomics integration reveal novel loci and biological insights into frailty,” *Nature Aging* (2025), DOI [10.1038/s43587-025-00925-y](https://doi.org/10.1038/s43587-025-00925-y).

**Subject:** Exact FinnGen R12 file or endpoint for the published HFRS GWAS

Dear FinnGen team / authors,

I am checking the provenance and availability of the FinnGen discovery summary statistics from the 2025 HFRS GWAS. The paper describes a custom continuous Hospital Frailty Risk Score based on 109 weighted ICD-10 codes, analyzed in FinnGen Release 12 (N=500,737), with a separate HFRS-without-dementia sensitivity analysis. I have checked the paper and its supplement, and the exact FinnGen endpoint code or file location for these custom-score GWAS results is not identified in the materials available to me.

Could you please confirm:

1. Whether the full HFRS and HFRS-without-dementia FinnGen R12 summary statistics are available in the standard public R12 manifest/browser, and, if so, the exact endpoint code and file path;
2. If these were run or deposited as custom GWAS results, the exact Sandbox analysis title and phenotype name used, the FinnGen release, and whether the results are listed in the signed-in `userresults.finngen.fi` browser or the release-specific `sandbox_custom_gwas` green-library location; and
3. If they are not available through either route, the authorized route to request or obtain the study-specific summary statistics.

For each available result, please provide the full genome-wide variant file (not only thresholded lead-variant tables) and source metadata needed for reproducible reuse: genome assembly, ancestry, effect scale/encoding, analytic sample size, variant filters, allele definitions, file checksum, and applicable license or terms. I will not substitute a related diagnosis endpoint for the custom HFRS score.

Thank you for your help.

Best,
[Name / affiliation]

## Evidence checked before drafting

- The [primary HFRS GWAS article](https://doi.org/10.1038/s43587-025-00925-y) reports FinnGen R12 discovery using the continuous HFRS, defined from 109 weighted ICD-10 codes, and a separate sensitivity analysis without dementia.
- FinnGen’s [R12 data-download documentation](https://finngen.gitbook.io/documentation/data-download) describes its public GWAS summary-statistics manifest and states that download instructions follow submission of the official access form. This general release route does not identify the custom HFRS file.
- The [FinnGen Handbook custom GWAS guide](https://docs.finngen.fi/working-in-the-sandbox/which-tools-are-available/untitled/custom-gwas-tool) says custom results can be located by their author-assigned analysis title in `userresults.finngen.fi` (sign-in required) and stored under a release-specific `sandbox_custom_gwas/<phenotype_name>` path. The paper and public materials reviewed do not identify that title or phenotype name. The local source audit did not verify an HFRS endpoint code or direct file route. No form has been submitted and no access request has been made.

This draft does not assert that the HFRS results are unavailable. It requests the exact authoritative file or the authorized path to it. No message has been sent.

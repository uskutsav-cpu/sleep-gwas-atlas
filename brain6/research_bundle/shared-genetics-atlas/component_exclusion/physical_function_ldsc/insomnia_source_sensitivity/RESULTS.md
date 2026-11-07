# Physical-function correlations: two insomnia sources on shared markers

The questionnaire-based Brain6 insomnia input and FinnGen R13 registry insomnia were each compared with UK Biobank walking pace and left grip strength using **1,020,637 common markers**. All four results use the same reference and unconstrained LDSC intercepts. This is a post-hoc source sensitivity, not independent replication or a matched removal of sleep items from FI or FFS.

| Physical measure | Insomnia source | r_g (SE) | p | Bonferroni p across these four correlations |
| --- | --- | ---: | ---: | ---: |
| Walking pace | Questionnaire | −0.2900 (0.0264) | 5.071×10⁻²⁸ | 2.028×10⁻²⁷ |
| Walking pace | Registry | −0.0608 (0.0319) | 0.0566 | 0.2264 |
| Left grip | Questionnaire | −0.1680 (0.0251) | 2.2765×10⁻¹¹ | 9.106×10⁻¹¹ |
| Left grip | Registry | −0.0947 (0.0345) | 0.0060 | 0.0240 |

The physical outcomes' source coding and allele conventions have not been independently reconstructed from the raw phenotypes. Assuming their usual increasing coding and the original insomnia coding, negative r_g means genetic propensity for higher walking pace or grip tends to correlate with lower insomnia propensity. These are genome-wide associations, not within-person correlations, causal effects or a specific shared locus. Individual p-values above are the rounded values printed by LDSC, multiplied by four for this limited family; this is not an atlas-wide correction.

To compare insomnia sources for the same physical outcome, a paired 200-block jackknife used the common markers. The point differences below use r_g values rounded in the LDSC log; confidence intervals are unadjusted.

| Physical measure | Questionnaire r_g minus registry r_g | Paired SE | Unadjusted 95% interval | Bonferroni p for two contrasts |
| --- | ---: | ---: | ---: | ---: |
| Walking pace | −0.2292 | 0.0374 | −0.3026 to −0.1558 | 1.85×10⁻⁹ |
| Left grip | −0.0733 | 0.0328 | −0.1376 to −0.0090 | 0.0508 |

The walking contrast crosses the two-contrast threshold. The grip contrast does not, and lies very close to the cutoff; its approximate calculation uses rounded point estimates. The result does **not** establish that questionnaire design alone caused either difference. The inputs differ in cohorts, ascertainment, ancestries, marker histories and model details. Questionnaire insomnia and the physical measures share UK Biobank participants; exact overlap and relatedness are unverified. The cross-trait intercepts for questionnaire insomnia are −0.1009 (SE 0.0070) for walking and −0.0510 (0.0066) for grip; they warrant attention to sample dependence and model fit but do not estimate participant overlap. FinnGen registry comparisons use a general-European LD reference, whose fit to Finnish LD is uncertain. Registry insomnia N is a study-level count rather than independently verified per-marker effective N; the questionnaire input uses an effective N. Original insomnia Z values were rounded in the archived source.

All source files were streamed to completion and checked for duplicate IDs; 1,020,637 common rows passed the stated non-palindromic allele compatibility rule and were retained by each LDSC fit. The shared physical-trait heritability delete vectors match exactly across the two insomnia fits, and the corresponding insomnia vectors match across physical traits. Each reconstructed jackknife SE matches its logged value to printed precision. This verifies the paired calculation and block alignment, not upstream phenotype coding, LD appropriateness or scientific validity. Common support is an already filtered subset, not all genome-wide variants. Logs include the existing compression runtime warning; the source gzip streams were separately read to completion.

The [pre-run protocol](PROTOCOL.md), [matching code](prepare.py), [input receipt](prepare_receipt.json), [run script](run.sh), [calculation](summarize.py), [results](results.json) and [delete vectors](rg_deletes.tsv) retain the exact procedure and hashes; complete per-pair logs are in this directory. Matched inputs are stored on the SSD. The earlier [registry-only pilot](../RESULTS.md) used a slightly larger support and a two-test family, so its numbers are not interchangeable with this four-input analysis. A valid comparison of physical frailty after removing exhaustion still requires a matched modified-score GWAS.

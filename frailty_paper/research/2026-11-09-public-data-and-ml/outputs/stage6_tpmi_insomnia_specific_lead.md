# TPMI insomnia-specific validation lead

On 28 October 2026, inspection of the [TPMI phenotype table](https://pheweb.ibms.sinica.edu.tw/phenotypes), filtered by phecode, identified an insomnia-specific phenotype that the earlier broad sleep-disorder screen missed. This corrects the impression that the available TPMI outcome is limited to broad sleep disorders. The earlier 12-variant inventory remains a screen of phecode 327, not insomnia.

| Phenotype | Cases | Controls | Case ICD codes shown by TPMI |
|---|---:|---:|---|
| [Insomnia, 327.4](https://pheweb.ibms.sinica.edu.tw/pheno/327.4) | 30,581 | 287,785 | 307.41, 307.42, 327.0, 780.52, F51.01, F51.02, F51.03, F51.04, F51.05, F51.09, G47.0 |
| [Organic or persistent insomnia, 327.41](https://pheweb.ibms.sinica.edu.tw/pheno/327.41) | 16,010 | 287,785 | 307.42, 327.0, F51.01, F51.03, F51.09 |
| [Sleep disorders, 327](https://pheweb.ibms.sinica.edu.tw/pheno/327) | 45,836 | 287,785 | 307.4, 327, 333.94, 780.5, F51.01, F51.02, F51.03, F51.04, F51.05, F51.09, F51.1, F51.3, F51.4, F51.5, F51.8, F51.9, G25.81, G47 |

The table shows the same control-exclusion code list for these three phenotypes: 307.4, 327, 333.94, 780.5, F51.01, F51.02, F51.03, F51.04, F51.05, F51.09, F51.1, F51.3, F51.4, F51.5, F51.8, F51.9, G25.81, G47. Subtypes can overlap; do not add case counts. Hospital-coded insomnia is a closer phenotype than the broader category but still differs in ascertainment and measurement from self-reported insomnia. [Chen et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12675292/) describe the hospital-based design and note possible undercapture of mild and common illness.

The [candidate inventory](stage6_tpmi_insomnia_candidate_inventory.csv) now records rounded results for both phenotypes for 12 of the 15 prototype markers. Three markers remain unverified, as in the earlier broad-sleep screen. The two phenotypes overlap and results are correlated. The generic site's effect-allele convention has supporting but not definitive production-source verification, so signs cannot yet be assigned to the prototype's counted allele. Two illustrative results from the inventory are shown below.

A [published-table comparison](stage6_tpmi_published_allele_comparison.md) subsequently found source-specific evidence for ALT-effect coding in selected examples, including a 327.41 insomnia record. These examples are not project markers and do not settle every source conversion or multiallelic contrast.

| Variant (GRCh38, displayed alleles) | Insomnia 327.4, displayed effect (SE), p | Organic/persistent 327.41, displayed effect (SE), p |
|---|---|---|
| [rs8067974, 17:45155199 A/G](https://pheweb.ibms.sinica.edu.tw/variant/17-45155199-A-G) | 0.0063 (0.0094), 0.50 | −0.013 (0.012), 0.29 |
| [rs11790060, 9:93440650 T/C](https://pheweb.ibms.sinica.edu.tw/variant/9-93440650-T-C) | −0.019 (0.011), 0.040 | −0.017 (0.012), 0.18 |

Among the 12 observed markers, rs11790060 is the only nominal p < 0.05 for insomnia 327.4 (p=0.040); none of the 12 is nominally below 0.05 for 327.41. These screens are postselected, unadjusted, and not evidence of score replication. Non-significance is not evidence of no effect. Rounded browser outputs are not a substitute for full precision, verified alleles, target LD, or a prespecified full-variant analysis. Participant/relatedness overlap with all discovery sources is not established. Neither subtype measures frailty or establishes score prediction in individuals. Next, obtain the official insomnia summary statistics, resolve the remaining three variants and source allele orientation, and evaluate what the summary data can support; full individual-level score validation still requires linked genotypes and relevant outcomes. The separate 22-marker score has not been fully assessed here.

TPMI requests acknowledgement and citation of [Chen et al.](https://doi.org/10.1038/s41586-025-09350-y) and [Yang et al.](https://doi.org/10.1038/s41586-025-09680-x).

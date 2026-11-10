# TPMI sleep-disorders lookup: rs7734045 and rs1032524

Two additional 15-marker prototype variants were found in the public TPMI PheWeb. Source positions use GRCh38 and match positive-strand Ensembl liftover of the prototype's GRCh37 coordinates. The allele pairs match, but TPMI effect-allele coding is not independently established here, so no signed comparison with the prototype weights is made.

| Variant | Project counted allele | TPMI GRCh38 position and allele pair | Browser effect (SE) | P |
|---|---|---|---:|---:|
| [rs7734045](https://pheweb.ibms.sinica.edu.tw/variant/5:136202653-A-C) | C | chr5:136202653 A/C | 0.015 (0.0079) | 0.062 |
| [rs1032524](https://pheweb.ibms.sinica.edu.tw/variant/1:201821698-T-C) | T | chr1:201821698 T/C | 0.016 (0.0096) | 0.090 |

The [sleep-disorders phenotype](https://pheweb.ibms.sinica.edu.tw/pheno/327) has 45,836 cases and 287,785 controls. These browser values are rounded; the broad diagnostic group is not insomnia alone. Neither result is statistically clear at nominal 0.05, and neither establishes a null effect. They were checked after candidate selection, alongside other markers, without a prespecified multiplicity adjustment. TPMI and all discovery/selection datasets have not been formally checked for participant or relatedness overlap. Different ancestry and phenotype definitions limit comparability.

These are exploratory single-marker results, not validation of an aggregate score or its discrimination and calibration in individuals. The browser labels SMIM32 and NAV1 as nearest genes, which does not establish that either mediates the associations. Previous lookups are linked in `stage6_tpmi_sleep_source_lead.md`.

TPMI requests acknowledgement and citation of [Chen et al. (2025)](https://doi.org/10.1038/s41586-025-09350-y) and [Yang et al. (2025)](https://doi.org/10.1038/s41586-025-09680-x).

## Additional marker: rs4886140

The [TPMI variant page](https://pheweb.ibms.sinica.edu.tw/variant/13:59259385-A-G) lists rs4886140 at GRCh38 chr13:59259385 A/G. Its project record counts A at GRCh37 chr13:59833519; the existing liftover audit matches the GRCh38 position on the positive strand. The rsID is multi-allelic in another reference, and TPMI effect-allele orientation remains unconfirmed. The browser displays effect 0.022 (SE 0.011), p = 0.052 for the same broad sleep-disorders phenotype and sample size. These rounded, postselected values do not show a statistically clear association at nominal 0.05 or establish a null effect, and they do not validate the score. DIAPH3 is the site's positional label, not established mediation. The same overlap, ancestry, phenotype and multiple-testing limitations apply.

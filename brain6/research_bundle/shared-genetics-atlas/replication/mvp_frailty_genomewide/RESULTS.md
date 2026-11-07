# MVP insomnia and three frailty measures: genome-wide sensitivity

Using the MVP European-ancestry insomnia source and the existing frailty inputs, exploratory LDSC estimated the following genome-wide genetic correlations on 1,103,239 valid SNPs in each comparison:

| Frailty measure | Genetic correlation | SE | p-value | Observed-scale h² (SE) | h² intercept (SE) | Cross-trait intercept (SE) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Frailty index (FI) | 0.5531 | 0.0390 | 1.0587×10⁻⁴⁵ | 0.1088 (0.0051) | 1.0235 (0.0091) | −0.0018 (0.0070) |
| Fried frailty score (FFS) | 0.3997 | 0.0352 | 8.205×10⁻³⁰ | 0.0623 (0.0027) | 1.0265 (0.0105) | 0.0180 (0.0076) |
| Hospital frailty risk score (HFRS) | 0.5697 | 0.0314 | 1.5712×10⁻⁷³ | 0.0542 (0.0024) | 1.1882 (0.0141) | 0.0013 (0.0074) |

The MVP observed-scale h² was 0.0489 (SE 0.0033), with intercept 1.1366 (0.0093). All three p-values are below the exploratory three-comparison Bonferroni threshold 0.0167; this does not correct for the wider project search. Elevated MVP and HFRS univariate intercepts require caution and cannot be assigned a cause from this analysis. The general-European LD reference may fit Finnish LD imperfectly. These estimates do not imply causality, shared causal variants, or clinical prediction. The estimates have not been formally compared against one another; their errors can be correlated, and differences in population and measurement also matter. Participant and relative overlap is unknown.

The FI includes insomnia and tiredness, FFS includes exhaustion and HFRS includes a somnolence-related diagnostic category. Consequently, this is not a test of physical frailty independent of sleep or exhaustion. Nor is it established independent replication. Compare with the prior [questionnaire-insomnia analysis](../../global_audit/brain6_frailty_covariance/three_frailty_common_support/RESULTS.md) only descriptively: source ascertainment and marker support differ.

The three original frailty inputs had identical SNP order and allele orientation across 1,134,953 rows, and their compressed hashes matched their earlier support record. LDSC merged 1,103,239 valid SNPs for every pair. The MVP stream was anchored to the archived ADHD/HapMap support; it read the complete source and its [completion receipt](stream-20261029T052003-99497.json) exactly matched the earlier audited receipt. Source p-values and rounded odds-ratio directions were used to reconstruct Z; variants with rounded OR=1 and duplicate eligible rsIDs were excluded, and the working constant effective N is approximate. These and other caveats are detailed in the [MVP input audit](../mvp_finngen_genomewide/RESULTS.md). LDSC issued a compression warning for the archived frailty files; hashes and successful read are useful checks but do not establish original processing provenance or model calibration.

See the [protocol](PROTOCOL.md), [run script](run.sh), [input hashes](input-hashes-20261029T052003-99497.txt) and [full LDSC log](ldsc-20261029T052003-99497.log). The computation was completed before the results were interpreted. No AlphaGenome output was used to train a model.

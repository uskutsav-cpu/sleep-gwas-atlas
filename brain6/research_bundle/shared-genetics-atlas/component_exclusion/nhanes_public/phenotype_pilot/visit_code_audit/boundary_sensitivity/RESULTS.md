# Sensitivity to the 4–5-visit category

The 2013–14 NHANES response category **4–5 healthcare visits** straddles the working FI cutoff. Assigning all of that category the lower score (0.5) or higher score (1) changed the FI for **263 of the same 2,433 primary respondents**. Neither assignment reconstructs their actual visit counts.

| Comparison | Lower assignment r | Higher assignment r | Higher minus lower (95% survey interval) |
| --- | ---: | ---: | ---: |
| Full FI vs adapted Fried score | 0.56416 | 0.56282 | −0.00134 (−0.00226 to −0.00042) |
| Full FI vs sleep symptom | 0.31078 | 0.31155 | 0.00077 (−0.00047 to 0.00200) |
| FI without sleep and fatigue vs sleep symptom | 0.23804 | 0.23904 | 0.00099 (−0.00031 to 0.00230) |

The assignments are two extreme working choices for that category, not confidence bounds on the true correlation or a validated scoring rule. Both assign a nonmissing value, so they alter the FI value but not its denominator for that item; the adapted score and eligible respondents remain fixed. The estimated shifts are small in this sample; that does not resolve the authors' original coding, changed question wording between survey cycles, other missing-item assumptions or selection. The intervals are post-hoc paired survey intervals, without correction for the broader atlas. These are observational phenotype correlations and do not provide item-excluded GWAS or evidence of a causal mechanism.

The [protocol](PROTOCOL.md), [comparison code](compare.R), [higher-assignment scoring copy](score_high.R), [counts](counts.tsv), [results](comparison.tsv), and [log](run.log) preserve the calculation. All 32 source hashes matched. Eligibility was checked in memory, score changes were restricted to the relevant category, and both correlations were checked against direct weighted calculations. No person-level records were saved.

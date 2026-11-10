# Shared association follow-up: insomnia and generalized weakness

The publicly available All by All results previously queried include exact position/allele rows for all 33 distinct project variants in both insomnia (`NS_333.2`) and generalized weakness (`SS_807.21`), separately in EUR and AFR genetic-similarity groups. We asked a narrower question than whether the scores predict disease: **is there association evidence for both coded outcomes at the same candidate?**

For each candidate and group, we used an intersection–union p-value `max(p_insomnia, p_weakness)`. The null is that at least one of the two trait associations is null. If both source p-values are valid, this max-p test is valid without assuming the outcomes or samples are independent. It does not test whether effects point in the same direction, whether they share a causal variant, or whether a score predicts either outcome. We additionally applied BH across the 33 candidates separately in each group; its standard FDR guarantee depends on the joint dependence structure, and it does not account for the much wider project search.

| Genetic-similarity group | Candidates | Both source p-values < 0.05 | BH q < 0.05 within 33 |
|---|---:|---:|---:|
| EUR | 33 | 0 | 0 |
| AFR | 33 | 0 | 0 |

The smallest max-p was 0.0581 in EUR and 0.282 in AFR. The lack of significant tests is not evidence that no shared effects exist; the groups have different sample sizes and LD, and weakness has a limited case count. No between-group test was performed. The candidates and outcomes were postselected from a large body of prior analyses, and discovery/target participant and family overlap remain unverified. These coded conditions are neither FI nor Fried frailty; their presence in the same biobank does not prove the same participants have both conditions or that the phenotypes are independent. Associations with both outcomes, had any been detected, would still require allele-direction and colocalization analyses before suggesting a shared causal mechanism.

The [66 candidate/group rows](stage6_aou_insomnia_weakness_conjunction.csv) preserve each original p-value, max-p and within-group BH q-value. Reproduce with `work/stage4_env/bin/python work/aou_complete/insomnia_weakness_conjunction.py`. The [receipt](stage6_aou_insomnia_weakness_conjunction_receipt.json) records source and output hashes, and the script checks unique candidates and identical position/REF/ALT pairs across the two source rows. These checks do not establish beta orientation, target-sample LD or score validation. Source and phenotype definitions are documented in the [AoU retrieval review](stage6_aou_full_candidate_association_review.md).

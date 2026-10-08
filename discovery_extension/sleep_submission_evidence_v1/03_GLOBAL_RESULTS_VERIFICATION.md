# Global numerical verification

Audit date: 2026-10-08. This technical report documents executed arithmetic checks, not a manuscript section.

The 100 external identifiers and 12 sleep identifiers define an exact Cartesian family of 1,200 distinct pairs. Both the historical matrix and pair-universe table match it with no omitted or duplicate pair. All 100 historical extension h2 records satisfy the recorded h2 Z≥4 and intercept≤1.2 gates. Historical files are read-only inputs to this audit. The original standalone sleep-trait h2/QC artifact remains missing. The eight retained replication rg logs provide first-pair sleep h2 records; all eight printed h2 Z/intercept values pass the recorded gate, but they are historical validation-log diagnostics and do not replace the original 12-trait discovery QC. All 13 standalone replication-source h2 logs also agree with their table; 9 sources pass and 4 fail the recorded gate.

The standard-library implementation recalculated BH from all 1,200 original P values. Every adjusted value agrees with the historical table; 603 comparisons remain significant at FDR<0.05. Separately replacing P by the two-sided normal tail of printed Z, or of printed rg/SE, yields the same 603 pair labels. Assuming nearest four-decimal printing, conservative joint rg/SE rounding bounds also yield exactly the same 603 labels. These are precision robustness checks; the original BH decisions remain authoritative. A retrospective 1,200-pair Bonferroni threshold gives 337 significant comparisons and does not replace the prespecified BH analysis.

The historical display is not an unrounded scientific result: rg, SE and Z are generally printed to four decimals, and the original discovery P field combines LDSC scalar precision conventions. `p_from_printed_z` and `p_from_printed_rg_se` are full floating-point calculations conditional on rounded input, not recovered full-precision P values. Normal 95% intervals are calculated from printed estimates and are not clipped at ±1.

All 217 replication identities agree with the two candidate/manifest locks and their original discovery estimates crosswalk to the global matrix. The exact threshold is 0.05/217 = 0.0002304147465437788. Historical displayed alpha is rounded to 0.000230414746544. Forty-one pairs were estimated, 17 were excluded by h2/intercept gates, and 159 had no qualifying source in the recorded search. All 41 estimates, SE and displayed Z values agree with the retained replication logs; normal tails reproduce their historical P values. No direction failures occur. Twenty-three cross the frozen threshold and 18 remain concordant but nonsignificant. Printed rg/SE and conservative rounding bounds do not change any of these 41 decisions.

All 23 successes reuse the discovery sleep GWAS and are labeled `EXTERNAL_OUTCOME_SIDE_REPLICATION`. They cover eight external phenotypes and seven sleep traits: they do not represent 23 independent diseases, loci, or fully independent two-trait replications. Eighteen use exact external definitions and five use `COMPARABLE_WITH_DOCUMENTED_DIFFERENCES`. Further source qualifications are in report 04.

A second implementation computes normal tails through the regularized incomplete gamma function, using a convergent series or continued fraction instead of erfc, and computes BH through explicit suffix minima. It checks all 1,200 rows and all 41 estimated replication rows; relative tolerance is 3×10⁻¹², absolute tolerance 10⁻³⁰⁰. The 20 source-free tests include missingness, duplicates, denominator shrinkage, extreme tails, invalid P/SE, covariance bounds, and unconstrained intervals.

| Gate | Status | Boundary |
|---|---|---|
| Complete 1,200/217 universes | PASS | Locked identifiers match |
| Original-P BH arithmetic | PASS | All 1,200 values checked |
| Replication threshold and classification | PASS_WITH_QUALIFICATION | All23 outcome-side only |
| Retained replication printed-log concordance | PASS | 41 estimates, no native replay |
| Independent numerical implementation | PASS | Source-free arithmetic |
| Original standalone sleep h2/QC for all 12 traits | BLOCKED | Missing original artifact; 8 validation-log records are partial evidence |
| Native discovery/replication LDSC rerun | BLOCKED | No eligible complete source/reference environment in this checkout |
| Unrounded jackknife estimates recovered | BLOCKED | Printed output cannot recover omitted precision |

Reproduce from repository root:

```sh
python3 discovery_extension/sleep_submission_evidence_v1/scripts/verify_statistics.py
python3 discovery_extension/sleep_submission_evidence_v1/scripts/independent_statistics.py
python3 -m unittest discover -s discovery_extension/sleep_submission_evidence_v1/tests -p test_statistics.py -v
```

Outputs: `tables/global_1200.tsv`, `significant_603.tsv`, `replication_family_217.tsv`, `replicated_23.tsv`, `numerical_summary.json`; receipts in `logs/numerical_execution.json` and `numerical_independent.json`. Hash verification of historical artifacts is a separate provenance gate in reports 01–02, and never substitutes for native reproduction.

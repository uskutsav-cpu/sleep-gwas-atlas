# Negative outcomes and blocked analyses

Audit date: 2026-10-08. Absence of an estimated result remains missing, never numerical zero.

All1,200 discovery comparisons are retained in `tables/global_1200.tsv`. The597 that do not meet BH-FDR<0.05 are negative discovery outcomes, not confirmed biological nulls. The selected panel has no documented negative-control design establishing that particular pairs should be zero. No synthetic or presumed-null phenotype was introduced. The original local queue contains 814 distinct rows: 217 selected priority candidates and 597 `PENDING_GLOBAL_NULL_CURATION` rows. The audit asserts exact identity equality between those 597 rows and the 1200-minus603 nonsignificant discovery set, and between the 217 selected rows and the complete replication family. These are post-global-result secondary candidates, not prespecified negative controls. The queue remains a source-free prioritization receipt, not an executed local experiment.

Within the complete217-candidate replication family,41 estimates were obtained. Twenty-three cross the frozen threshold;18 are directionally concordant but nonsignificant; no estimated pair has direction disagreement. Seventeen candidates had sources failing h2/intercept eligibility, and 159 had no qualifying dataset in the documented pre-result search. The latter labels describe the historical search boundary; they do not establish that no suitable study exists anywhere today. In the new verification table they retain the historical labels and exact search fields. The17 QC-stopped and 159 source-unavailable rows have `NA` estimates and P values. They are not combined with the18 tested nonsignificant outcomes as failed biological replications.

No equivalence margin was prespecified, and no power calculation tied to a scientifically justified smallest effect was available. A nonsignificant P or overlapping confidence interval is not evidence of equivalence. Finite positive correlation estimates outside a rejection threshold cannot be translated into absent shared biology.

| Component | Status | Required restart input |
|---|---|---|
| Source-free global/BH/replication arithmetic | PASS | Retained processed tables |
| Printed replication log concordance | PASS | Retained41 printed results |
| Full-precision original estimates | BLOCKED | Unrounded outputs or eligible native replay |
| Native discovery and replication LDSC replay | BLOCKED | Exact eligible sleep/outcome summary statistics, ancestry-matched LDSC references and reproducible original runtime |
| Calibrated effect-difference inference | BLOCKED | Joint block-jackknife estimates or defensible sampling covariance |
| Smoking/adiposity conditional analysis | BLOCKED | Eligible multivariate genetic/sampling covariance, frozen model and source eligibility |
| Formal correlated-phenotype cluster inference | BLOCKED | External phenotype genetic and sampling covariance matrices |
| Prespecified negative controls/equivalence | NOT_AVAILABLE | Defensible prespecified design and relevant precision/power |
| Local correlation / PLACO+ / fine-mapping / colocalization | BLOCKED_UPSTREAM | Verified dense variant statistics, ancestry-matched LD and method-specific QC |
| Causal, gene, therapeutic or molecular mechanism claims | NOT_SUPPORTED | Appropriate additional design and evidence |

Blocked local or mechanism analyses are encoded as unavailable, not zero loci or zero shared variants. HapMap3-munged LDSC inputs are not locus-complete dense data. The source-free robustness calculations executed here add no native GWAS experiments, causal effects or new loci. Future source acquisition and joint-covariance work must use a versioned namespace and retain all failures without replacing the locked historical outputs.

# Genome-wide correlation recheck

## Insomnia–ADHD

The 396-row source atlas includes exactly one insomnia–ADHD row: r_g=0.4033, SE=0.0318, p=5.836×10⁻³⁷, BH q=2.56784×10⁻³⁵. Reapplying Benjamini–Hochberg to **all 396 stored p-values** reproduces all stored q-values exactly, including this result (ninth-smallest p). This is an independent check of the table's multiple-testing arithmetic. A later [same-data LDSC re-execution](insomnia_adhd_rerun/RESULTS.md) matches the printed r_g, SE and p using archived inputs that match the later Brain6 freeze. The original atlas log and complete original-run provenance remain unavailable.

Reported r_g and SE are rounded. Their ratio gives z=12.6824 and a normal p≈7.40×10⁻³⁷, whereas the atlas stores z=12.7 and p=5.836×10⁻³⁷. The discrepancy is compatible with precision lost in reporting but cannot be fully resolved without the original unrounded run. Do not recompute the published BH correction from the rounded r_g and SE or imply that the p-values match exactly.

## FI–FFS

The archived native FI–FFS LDSC record gives r_g=0.7534110683, SE=0.0199216376 and 1,157,249 regression SNPs. I independently recalculated the jackknife SE from its 200 delete estimates using sqrt[(199/200) Σ(delete_i − mean(delete))²], obtaining 0.01992163760340901 versus the archived 0.019921637603409075. That calculation verifies the stored uncertainty arithmetic. A later [same-data standard LDSC re-execution](fi_ffs_rerun/RESULTS.md) also reproduces the reported r_g, SE, z and SNP count; it does not establish genetic causal inference, source allele convention or an independent cohort. The archived p is zero due to numerical underflow; it is not literally zero. Exact sample overlap within UKB is unknown.

`atlas_396_source.tsv` is the project’s original `results/atlas/trait_pairs.tsv`. `fi_ffs_native_exact.json` is copied from `/Users/swethasunilkumar/Documents/Codex/2026-10-08/g/work/sleep-gwas-atlas/frailty_paper/full_translational_psychiatry_research_v3/logs/native_rg_fi__ffs.exact.json`. `check_global.py` reproduces the audit and records source-copy hashes in `audit_receipt.json`.

# Completed PLACO+ scan audit

This audit is prepared while the scans are running. `audit_completed.py` refuses to read `.partial` files and requires both final pair files. It traverses each final result together with the two original aligned input files and checks exact row ordering, SNP and effect-allele identity, rounded Z values, marginal normal p-values, uniqueness, status fields, finite valid PLACO+ p-values and full support counts. It requires all 2,191 precomputed control SNPs to match their original allele-aligned Z values and independently computed PLACO+ outputs in each pair. Numerical failures and excluded extremes remain in the original 2 × 1,138,262 correction family. It hashes all input and completed output files, tabulates statuses, and saves every marker–pair crossing 0.05/(2×1,138,262), alongside both marginal p-values and the max-p conjunction result. It does not claim distinct loci, apply an LD-clumping method, validate the nuisance model, or resolve source model or precision concerns. The threshold-crossing file is diagnostic and exploratory, not a discovery list.

After both terminal jobs finish successfully, run from the workspace root:

```sh
python3 outputs/shared-genetics-atlas/placo_diagnostic/audit_completed.py \
  --input '/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support' \
  --results outputs/shared-genetics-atlas/placo_diagnostic \
  --out outputs/shared-genetics-atlas/placo_diagnostic/completed_audit \
  --controls outputs/shared-genetics-atlas/placo_diagnostic/control_results.tsv
```

`test_completed_audit.py` verifies a small synthetic two-pair run with a numerical failure, the full-family threshold, a mismatched SNP and a missing completed file. The fixture validates the audit logic, not PLACO+ calibration. Successful execution of the audit itself still does not prove valid GWAS marginal tests, proper population-matched LD, independence of samples, shared causality, or independent replication. Check source-specific full-precision statistics for every threshold crossing and compare candidate rankings under alternative justified input models.

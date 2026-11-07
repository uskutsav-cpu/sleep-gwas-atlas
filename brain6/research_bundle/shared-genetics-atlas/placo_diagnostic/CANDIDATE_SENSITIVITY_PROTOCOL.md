# Follow-up of completed PLACO+ threshold crossings

Prepared before the full scan finished. After running the completed-output audit, join its threshold-crossing SNP–trait rows to the original filtered R13 insomnia, raw Atkins FI and Ye FFS association sources by unique rsID and effect/other alleles. A direct or reversed allele pair is accepted and the beta is aligned; missing, mismatched, duplicate valid rows, invalid values, zero beta with non-unit p, and a source-p-derived Z that disagrees with the stored three-decimal Z by more than 0.000501 are held rather than silently treated as zero. The candidate file must match the SHA-256 in `completed_audit.json`. Save the source beta, SE, p, orientation, reconstructed signed source-p Z, Wald Z and status for each pair. The FI source beta/SE are rounded and the FFS effect-allele adapter and exact BOLT model mapping remain unconfirmed. A unique rsID and allele pair alone does not verify genome build or all multiallelic representations.

For eligible rows, evaluate official PLACO+ with source-p-derived Z and beta/SE-derived Z with the same previously fitted full-support nuisance parameters and numeric settings. Values are sensitivities, not independent or fully recalibrated tests. A model-appropriate result requires full precision, source-method review and nuisance parameter re-estimation; do not treat the Wald result as automatically better. Select no new candidates using these sensitivities. Even stable threshold crossings need ancestry-matched LD, shared-causal evidence and independent two-trait replication. Empty, failed or held inputs are not null evidence.

Run only after `audit_completed.py` succeeds:

```sh
python3 outputs/shared-genetics-atlas/placo_diagnostic/candidate_source_audit.py \
 --audit outputs/shared-genetics-atlas/placo_diagnostic/completed_audit \
 --s '/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/filtered/R13_INSOMNIA_HM3.tsv' \
 --fi '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/data/raw/frailty.txt.gz' \
 --ffs '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv' \
 --out outputs/shared-genetics-atlas/placo_diagnostic/completed_audit
Rscript outputs/shared-genetics-atlas/placo_diagnostic/candidate_sensitivity.R outputs/shared-genetics-atlas/placo_diagnostic/completed_audit
```

`test_candidate_source_audit.py` checks allele reversal, exclusion of duplicate and missing matches, R sensitivity evaluation, and candidate file checksum enforcement on a synthetic fixture. These tests verify workflow behavior, not GWAS calibration or biological validity.

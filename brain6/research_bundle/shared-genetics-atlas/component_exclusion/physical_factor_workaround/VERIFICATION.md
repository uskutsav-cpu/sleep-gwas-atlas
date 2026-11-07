# Verification notes

## Reconnected-SSD rerun, 6 November 2026

The SSD became available again. A [ten-file hash audit](ssd_reconnection_hash_audit_20261106.json) confirms that the five retained source inputs and five common-support inputs match the archived receipts. The [retained-input LDSC replay](replay_ldsc.py) then reran all nine physical-indicator and insomnia comparisons in a new SSD directory. All nine printed genetic correlations and SEs matched the original logs; all 27 vectors of 200 heritability/covariance block-delete values matched with maximum absolute difference zero. Each run retained 1,053,614 SNPs. Hashes of the five common-support inputs and 88 reference/code files were checked before and after the replay. The [receipt](ldsc_replay_20261106/receipt.json), nine logs and 27 delete files are copied into this bundle so the verification record remains available if the SSD is disconnected.

This extends the earlier saved-summary reproduction to the LDSC calculations from retained munged common-support files. It does not rerun original raw GWAS production or initial munging, establish phenotype provenance beyond the existing source audits, resolve overlap/reference limitations, supply independent validation, or reconstruct exhaustion-deleted frailty. The same inputs and reference are deliberately reused; matching results are reproducibility evidence.

## Earlier checks

The three physical-source SHA-256 hashes match the earlier archived receipts. The questionnaire source's SHA-256 was `07d47cd27d616885b4143b66551b252aca50cb0671b55f0a23ccdc7105ca9d13` in both the retained download and a separate complete second retrieval; no publisher checksum was identified in this check. The preparation scripts read gzip files through completion, checked ordered unique SNP IDs and valid values, and retained the same 1,053,614 compatible nonpalindromic SNPs. Each LDSC log reports that number after allele QC.

The evaluation reconstructs 200-block genetic-correlation standard errors from the printed delete values and checks them against the LDSC logs to within 0.00005. It also checks that each trait's heritability delete vector is identical across the paired fits. An independent linear solve of the recorded GLS covariance matrices reproduced both point estimates and weights to numerical precision (largest weight difference 1.8×10⁻¹⁸). These checks verify inputs and calculations, not participant overlap, source phenotype coding, factor-model adequacy, independent replication or clinical utility.

# Frailty local-sharing readiness audit

Audit date: 2026-09-24 (UTC)

Plan: `frailty_paper/config/analysis_plan_v1.yaml`

Conditional decision: D07 in `frailty_paper/config/analysis_plan_decisions.tsv`

## Decision

Do not start a frailty-specific LAVA association run yet. The technical heritability gate passes for the Frailty Index and all 12 locked sleep traits, the 13 full dense harmonized files are readable and match their QC row counts, and the pinned LAVA reference payloads pass direct SHA-256 verification. The project-level Phase 9 gate is still unmet: no prioritized sleep–frailty association has independent replication, exact participant intersections remain unknown, and D07's complete frailty-specific family has not been frozen. The existing LDSC munged files are HapMap3-filtered and do not meet the frozen dense-input rule for primary locus work.

If local sharing is later justified without exact participant intersections, it must be explicitly sensitivity-only, use overlap estimates rather than assume zero overlap, and preserve a predeclared complete family. This audit contains no LAVA association results.

## Subsequent locked sensitivity-family execution update (2026-09-24 UTC)

The earlier D07-era decision above is historical and is superseded for this narrowly defined sensitivity analysis by D13 and `frailty_paper/config/lava_frailty_sensitivity_v1.yaml`. D13 froze all 12 sleep traits × 2,495 regions, the overlap-estimation method, inherited local-univariate threshold, all-slot BH family, and completeness gates before frailty-specific LAVA result access. The unmet independent-replication gate still bars confirmatory/Phase 9 Tier 1 local claims; execution is permitted only as sensitivity analysis.

The complete 13-file source family was rechecked against its locked byte counts, SHA-256 values, schemas, row counts, and nonempty coverage across chromosomes 1–22. The chromosome-sharded preparation completed with 13 traits, 286 shards, 12 pair-specific overlap matrices, 264 chromosome input-info files, 2,495 loci, and 29,940 candidate bivariate slots. Its external `input_manifest.json` SHA-256 is `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`; the read-only preparation verifier passed. The 44 pinned reference payloads and exact R 4.3.3/LAVA 0.1.5 implementation validator passed.

The per-locus runner is in progress. A verified pilot receipt for actigraphy sleep duration × Frailty Index at locus 1 records FI local p=0.185375 and a negative local variance/drop for the sleep phenotype; this is a QC/status example, not a bivariate result or family-level inference. The runner preserves separate per-trait statuses so a dropped phenotype does not discard the other valid univariate result. At the time of this update, the live run had verified 15 receipts, with no launch failures in the active invocation; no final family tables or scientific conclusions are available yet. See scripts `43_prepare_lava_sensitivity_inputs.py`, `44_run_lava_sensitivity_locus.R`, `45_run_lava_sensitivity.py`, and `46_collate_lava_sensitivity.py`.

## Interim family-completeness gate check (2026-09-24 09:07 UTC)

The same live run had 381/29,940 receipts verified at this checkpoint, all for
the first pair (`accel_sleep_duration`). Among these, 94 receipts had
`process_status=PROCESS_FAILED` with `process.locus returned NULL`. Under the
collator's frozen locus-failure definition (a locus index is flagged if any
pair has a process failure or either per-trait univariate status other than
`TESTED`), 271 distinct locus indices were already flagged in this first pair.
That is 10.86% of the 2,495-locus family, above the locked maximum of 1% (at
most 24 loci). The family-completeness gate is therefore already unable to
pass under the current collator and unchanged lock, even if every remaining
receipt passes. This is an interim QC failure, not a completed family result
or evidence of absence of local sharing. The resumable run remains active so
all fixed slots and failure/ineligibility states can be collated; no local
result may be interpreted if final family QC fails. The threshold and failure
classification are not being changed after result access.

## Receipt-failure diagnosis at 561 verified slots (2026-09-24 09:37 UTC)

At the runner state checkpoint (`last_trait=accel_sleep_duration`,
`last_locus_index=561`, 561/29,940 verified receipts, zero launch failures in
the active invocation), the frozen collator's first-pair locus rule flagged
394/2,495 loci (15.79%). This is a lower bound on family-level failure and is
already above the locked maximum of 24 loci (1%). The count decomposes into 135
loci where `process.locus` returned `NULL`, 164 processed loci where the sleep
phenotype was dropped, and 95 processed loci where frailty was dropped. These
categories are disjoint at this checkpoint.

A read-only check of all affected per-locus receipts and corresponding R logs
found the expected LAVA diagnostic in every case: the 135 `NULL` results report
negative variance for all phenotypes; the 164 sleep drops and 95 frailty drops
report negative variance for the respective phenotype. No unexplained worker
or reference-path failures were present in these canonical receipts. Earlier
reference-path launch attempts for locus 1 are retained in its log, but the
canonical receipt was subsequently processed successfully and is not counted
as a process failure. This diagnostic supports treating the high failure
fraction as the locked QC outcome; it does not justify changing thresholds or
interpreting local sharing. The complete fixed family remains in progress for
transparent collation.

Reproduction: `python frailty_paper/scripts/47_audit_lava_interim_receipts.py --through-index 561 --out frailty_paper/analysis/lava_interim_receipt_diagnostic_2026-09-24_561.json`. The read-only auditor verifies each receipt/log identity and diagnostic, applies the same failure predicate as the frozen collator, and records hashes for the receipt/log prefix and relevant scripts/config. The resulting JSON reports a receipt/log manifest SHA-256 of `675833d3115248be83414f65c171e3f9b9f629885f1ab89fdeabf3d1458ff154` and the unchanged lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.

## Heritability and harmonized-input checks

`frailty_paper/results/frailty_v1/ldsc/h2_master.tsv` records PASS for the primary Frailty Index (Z=21.86, LDSC intercept 1.020) and every locked sleep trait (sleep-trait Z range 9.69–28.53; intercept range 1.008–1.117). These satisfy the frozen inclusion gate (Z >= 4 and intercept <= 1.20). This is an audit of existing results, not a rerun.

For each of the 13 primary dense harmonized files on the designated external volume, a read-only pass verified the exact header `SNP CHR BP A1 A2 FRQ BETA SE P N`, full gzip decompression to EOF, row count against its `.qc.txt`, GRCh37/hg19 output build, and SHA-256. Existing QC sidecars record removal of ambiguous alleles, MAF <= 0.01, MHC variants and duplicate SNP IDs. Sleep apnea originated in GRCh38 and was lifted to hg19 with the recorded UCSC chain SHA-256 `14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1`.

| Trait | Rows | Compressed bytes | SHA-256 |
|---|---:|---:|---|
| Frailty Index | 6,291,635 | 106,450,123 | `241ed294ccc367b22787a0c051f8d377c020e16089e82920a555286f3148e6de` |
| Actigraphy sleep duration | 6,482,514 | 144,197,232 | `2df99672307d103ac00d6f4ee280ee4ed788c8c777abd39f041177b48ae43bc4` |
| Chronotype | 6,482,405 | 145,654,493 | `c561f7951762cc8e98b7969f220150fafbd11d00b0dd470f18e0a245adab2286` |
| Insomnia | 6,077,635 | 125,531,591 | `81c71d838b19f0297e09ec027d00382732fee1bfd112f44cfd2feca44a4f8bf5` |
| Long sleep | 6,549,769 | 147,859,672 | `150b4136404adb77a8484c6326a6f233290762a71ca4c851af9462de820d0f4c` |
| Napping | 7,042,618 | 155,409,529 | `d9346be15b3f84b2afef1a9a1ce754dc776b90cc3aa7e78dad18d63607714c43` |
| Short sleep | 6,549,791 | 146,959,493 | `47948b1e63eb15107545cde531e06f82227afce891f24b813039b930e574497d` |
| Sleep apnea | 6,841,149 | 167,251,776 | `a9c37682a1da7687277781518a2d4b94b19a7985270941e7209d594f5f353aa5` |
| Sleep efficiency | 6,482,540 | 144,204,520 | `18c766173bf0acca325b927ea4d7118ad51dd9e3310924132707dd12c0d34fb7` |
| Sleep timing | 6,482,540 | 144,192,547 | `69f3f10605aa07017d108ad8c1b5f9c02a95743ab507823e7bf4c9f50a77d2d4` |
| Sleep duration (continuous) | 6,549,809 | 145,984,615 | `f2c0b9eaa4623e9757de86d6542fb3267874a2def626ca548fccb3ab97fced41` |
| Daytime sleepiness | 6,549,823 | 144,902,667 | `1a7002ce3679ab15cd3a6af2ba7bd5d54b937e19693ee6b5b1ec1e7c75cd3afe` |
| Snoring | 7,168,629 | 163,291,287 | `b2a2e6871941e56c5eee3799b1af63c6933329b68a6d8def16d4b503c08c7683` |

Inputs reside under `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized/`. The current repository symlink is `frailty_paper/data`, which is unrelated to these analysis-workspace files.

## LAVA reference and existing tooling

The physical reference directory `/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1` contains 22 `.info` and 22 `.bcor` payloads. All 44 sizes and SHA-256 hashes match the pinned reference provenance (SHA-256 `35ab371935b9c260dab21be248f4907b693565787ef9f62334e7d8030204270e`) and extracted-manifest digest (`5e250575a5afba843d2bed12cbdc4cc87071267d3096ca1f4b25e50da9dba579`), total 15,653,917,776 bytes. The locus file has 2,495 regions and matches policy SHA-256 `462e81bce9ca85c9f11c0feb4d1bda24dd2c0fe095bfcc83fd32f77dff9a0882`.

Reproducible reference verification:

```sh
python frailty_paper/scripts/42_verify_lava_reference_mount.py \
  --reference-root '/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1'
```

It passed on 2026-09-24. The repository-wide `scripts/146_manage_lava_reference_archives.py --verify-state` still rejects `ref/lava/ukb_v1.1` because that repository path is a symlink; the new helper verifies the physical external directory against the already pinned provenance without changing that contract.

The established `scripts/02_munge.sh` uses `--merge-alleles .../w_hm3.snplist`; existing LDSC munged files therefore are not the full dense inputs required for local locus analysis. The full dense harmonized files above retain the required SNP, effect-allele, beta, standard-error and N columns. At the time of this initial readiness assessment, only the original 45-trait/396-pair LAVA setup existed and no frailty-specific runner or result family was available. The subsequent D13 sensitivity-family execution update above records the separately locked adapter and current run. The system `Rscript` checked during the initial audit could not load `GenomicSEM`; the dedicated run used the pinned R 4.3.3/LAVA 0.1.5 runtime instead.

## Replication and family gates

The current FI-versus-sleep table has 12 global LDSC estimates, with 9 retaining q < 0.05 under the inherited 396-pair family. That is global discovery evidence, not regional evidence. Phase 9 of the project brief prioritizes Tier 1 local analysis for replicated primary frailty relationships; Phase 8 and `analysis/sample_overlap_assessment.tsv` show expected UK Biobank sharing for 11 sleep sources and unknown exact participant intersections. The available LDSC cross-trait intercepts can inform overlap sensitivity analyses but do not turn these estimates into independent replication.

No relationship meets the Phase 9 Tier 1 priority. The original D07 conditional decision remains unmet for confirmatory local analysis. For sensitivity-only analysis, D13 superseded the D07-era status and froze all 12 sleep × 2,495 reference loci (29,940 slots), the overlap method, inherited thresholds and family gates before result access. The run continues unchanged for complete transparent collation.

At the read-only 2026-09-24 10:37 UTC checkpoint, the first pair had 993 contiguous receipts and 704/2,495 flagged loci (28.22%), exceeding the locked maximum of 24. This prefix includes 236 all-phenotype negative-variance failures, 279 sleep and 169 FI negative-variance drops, plus 20 process failures at loci 950–969. Each of those 20 logs reports that no specified SNP IDs are present in the reference. The 20 consecutive chr6 loci span 25,684,630–33,864,262; both prepared pair inputs contain zero variants in chr6:25,000,000–34,000,000, matching the pre-existing MHC exclusion recorded in all 13 harmonized-input QC sidecars. This makes the error consistent with the locked input exclusion leaving the reference loci empty; it is not a negative-variance outcome. The frozen collator still counts these process failures, and they are not removed from the denominator after result access. All 20 loci and exact diagnostics are retained in `analysis/lava_interim_receipt_diagnostic_2026-09-24_1038.json`.

The first-pair completeness gate has failed under the unchanged rule, so no family-level local-sharing inference is available. Continue the predeclared run for all 29,940 slots, preserve every status and error, and do not change thresholds, inputs, or locus eligibility post hoc.

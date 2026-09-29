# Brain6 confirmatory LAVA source-rescue boundary (2026-09-27)

This report records the prospective source-first rescue attempt against the
unchanged seven-trait, 2,495-locus family. It is updated only with separately
identified pilot receipts; pilot results cannot change a source-admission
decision. The pre-outcome burden was committed as `811f8332` before the new
source audits. The protected canonical result is 13,745 `TESTED`, 3,720
`NOT_RUN`, zero `FAILED`, with a frozen maximum of 873 `NOT_RUN`.

## Exact gate arithmetic

The family needs at least **2,847** additional tested cells. Perfect repair of
the exact long-sleep trait alone leaves **2,429** `NOT_RUN`. Even perfect repair
of every *other* trait leaves the **1,291** unresolved long-sleep cells, still
**418 above** the ceiling. The unique smallest set whose hypothetical perfect
repair can meet the gate is long sleep, insomnia, Parkinson disease, and MDD;
it leaves 533. These are upper bounds on possible recovery, not forecasts.
The deterministic source is
[`phase1/summary.json`](phase1/summary.json), with exact trait and pair
burdens in the adjoining TSVs. No bivariate LAVA result was generated while
the full-family gate failed.

## Trait-specific source boundary

| Trait | Frozen `NOT_RUN` | Source evidence and unresolved field or model |
|---|---:|---|
| Long sleep | 1,291 | Dashti 2019 ≥9 h versus 7–8 h, 34,184 cases and 305,742 controls. Nine local containers resolve to one 14,661,601-row, ten-field release with no `N`, `NMISS`, per-variant cases/controls, or export manifest. Broad and EBI are byte-identical. Primary and coauthored methods support BOLT-LMM as its analysis class, but the exact released statistic/model mapping and analyzed N are not established. Overall 339,926 is a documented LAVA proxy possibility, not a source-verified N; the frozen protocol admits it only for diagnostic/sensitivity use. The prior locked 88-locus primary linear pilot was 44/88 `TESTED`, equal to canonical at those loci, below the prespecified nine-net-locus advancement rule. A distinct binary arm cannot be chosen by its result. See [`public_longsleep/SOURCE_EVIDENCE_LEDGER.md`](public_longsleep/SOURCE_EVIDENCE_LEDGER.md), [`local_archive/LOCAL_ARCHIVE_AUDIT.md`](local_archive/LOCAL_ARCHIVE_AUDIT.md), and the prior pilot adjudication. |
| Insomnia | 671 | Jansen 2019 UKB-only release has literal per-SNP `N`, effect-allele `OR`, and GRCh37 positions. Historical harmonization discarded that N and wrote constant effective N. Official CNCR README is pinned in `other_traits/`. The exact regression command and per-SNP case fractions remain unavailable; a single predeclared native-N trait-only technical pilot may assess feasibility with the explicit cohort-fraction approximation, but cannot by itself admit a confirmatory bivariate input. |
| Parkinson disease | 636 | Nalls 2019 release has per-SNP `N_cases` and `N_controls` but mixes direct and UKB proxy cases across a heterogeneous meta-analysis. The raw file lacks rsID; source-specific count interpretation, allele-aware reference mapping, and pair overlap remain unresolved. |
| MDD | 589 | Howard 2019 public file has no SNP N, INFO, or coordinates. PGC MDD2025 has richer counts and a prior 108/111 numerical pilot, but its cohort-wise effective-N and case-fraction mapping remains on scientific hold; pilot success cannot select the model. |
| ADHD | 223 | Demontis 2023 provides variant-level case/control fields, but participating cohorts and case fraction vary by SNP. A source-supported meta-analysis-to-LAVA mapping is needed. |
| Bipolar disorder | 168 | Mullins 2021 PGC fields include effective rather than proven literal analyzed case/control counts; no total-N substitution is admitted. |
| Schizophrenia | 142 | Trubetskoy 2022 PGC effective-N conventions and case/control/trio structure require source- and method-supported mapping before any replacement. |

The [other-trait source ledger](other_traits/SOURCE_FIRST_LEDGER.md) identifies
the exact raw files, hashes, source fields, and direct-source versus inferred
semantics. Existing source and overlap manifests remain untouched.

The source-first insomnia v2 technical pilot was frozen and committed before
outcomes (`225e89ac`). It materialized 88 representative loci and 207,643
receipt-verified SNP rows with original SNP/alleles/Z and native source N,
then launched four disjoint 22-locus workers. Its immutable terminal receipt
records worker exits `[0,1,0,0]`; worker 2 failed the **predeclared zero-failure
condition**. The v2 pilot is therefore terminally failed and cannot authorize
an insomnia-only complete screen or any confirmatory promotion. A separately
versioned worker-2-only postmortem (`c5fb4d7c`) is diagnostic and cannot
change that decision. The postmortem reproduced a single failure at locus
**950**: its original and native-N input shards both had zero SNP rows, so
LAVA reported fewer than three shared SNPs. This is empty-input handling in
the pilot runner; it does not discredit the README's native-N semantics.
The three successful v2 worker files plus the separately labeled v3 diagnostic
file contain a **diagnostic-only** 88-locus tally of 47 `TESTED`, 40 `NOT_RUN`,
and one `FAILED`, versus 44/44/0 at the same canonical loci. The v2 terminal
failure is not retroactively reclassified or advanced. No source model,
threshold, or locus set was selected from these outcomes. The full explanation
is in
[`other_traits/INSOMNIA_NATIVE_N_WORKER2_POSTMORTEM_V3_RESULT.md`](other_traits/INSOMNIA_NATIVE_N_WORKER2_POSTMORTEM_V3_RESULT.md).

## Public, local, and controlled-access routes

The long-sleep search checked the original Broad/SDKP release, GWAS Catalog
GCST007560/EBI mirror, local SSD copies and historical archives, older Jones
and newer publications, public portal/consortium alternatives, and modern
phenotype-adjacent datasets. No second public full-genome source combines the
**same ≥9 h versus 7–8 h** contrast with auditable analyzed N and an exact
LAVA mapping. SleepChart (>8 h versus 6–8 h), ordered/continuous duration,
≥10 h, and PRS weights are different analyses. A 2025 UKB category analysis
is an acquisition lead, not an available full-GWAS file. The local audit found
no project receipt for an approved UK Biobank individual-level environment or
MVP/dbGaP access. That does not assert investigators lack institutional
access; it means no lawful, technically available access was documented here.

The [prepared narrow custodian request](../lava_longsleep_source_rescue_v1/LONG_SLEEP_DATA_REQUEST.md)
asks for variant-keyed analyzed N/NMISS or a source-specific proof of fixed N,
the exact released BOLT versus logistic export/model statement, and phenotype
and allele/QC details. It has not been sent. A complete exact-phenotype rerun
with these fields is another admissible acquisition path after authorization
and prospective method review. No guessed N is entered in a confirmatory
input.

## Current decision and work boundary

**TRUE EXTERNAL IMPOSSIBILITY with presently available source and access
evidence.** The frozen QC gate remains failed. No new source has been admitted
for a full-family confirmatory replacement. The exact-phenotype long-sleep
source lacks the critical N/model contract, all nine local containers and
authoritative public copies are the same payload, and no authorized
individual-level rerun environment is documented. Its 1,291 unresolved cells
exceed the 873 ceiling even if every other trait is repaired perfectly. Thus
the gate cannot be reached from the available admissible inputs. The
source-faithful insomnia pilot is technical evidence and cannot overcome that
floor. Therefore do not launch a full
seven-trait run, bivariate LAVA, or candidate promotion from the present
files. The original 25 pair-specific candidates and 20 reconciled regions
remain preserved and unpromoted. Confirmatory replication, GWAS fine-mapping,
trait–trait/GWAS–QTL colocalization, and enrichment **for promoted regions**
remain conditional. Their available bounded/exploratory feasibility and
evidence tables are complete in
[`../brain6_bounded_manuscript_v1/`](../brain6_bounded_manuscript_v1/).
The independent [Phase 23 bounded-package audit](bounded_final_audit/PHASE23_BOUNDED_PACKAGE_AUDIT_20260928.md)
verified 28/28 package hashes, 12/12 table input hashes, four/four method
input hashes, 25 candidate and 20 region joins, the canonical read-only
validators, and the four-worker execution lock. The current rescue-specific
[`INTEGRITY_AUDIT.json`](INTEGRITY_AUDIT.json) additionally checks frozen
anchors and the SSD pilot's immutable terminal/postmortem receipts.

## Exact restart after source information arrives

Use the existing
[`confirmatory_source_manifest.template.json`](../brain6_bounded_manuscript_v1/confirmatory_source_manifest.template.json)
and
[`signed_method_review.template.json`](../brain6_bounded_manuscript_v1/signed_method_review.template.json)
at new source-specific paths. Bind original source files, custodian/method
documents, exact phenotype, N/model mapping, build, alleles, case fraction or
continuous treatment, cohort overlap, hashes, scientific review, and an unused
SSD output root. The guarded command from the repository root is:

```sh
python3 brain6/scripts/brain6_confirmatory_workflow_v1.py --manifest /ABSOLUTE/PATH/TO/NEW_SOURCE_MANIFEST.json
```

That command must return ready before the separately reviewed, hash-bound
stage adapters can run with `--execute`. It fails closed on the current
placeholder manifest. The [restart protocol](../brain6_bounded_manuscript_v1/CONFIRMATORY_RESTART.md)
specifies the pilot, multi-trait projection, four-worker resumable family,
frozen gate, promotion, and downstream sequence. The old canonical, v2,
roundoff, pilot, PLACO, and bounded outputs must never be overwritten.

## Terminal outcome inventory

| Requested outcome | Current evidence-bound answer |
|---|---|
| Confirmatory LAVA family | 13,745 `TESTED`, 3,720 `NOT_RUN`, zero `FAILED`; frozen QC **failed**. No source-admitted replacement family exists. |
| Admitted repairs | None for confirmatory replacement. Jansen native per-SNP N is source-verified and was used only in the versioned sensitivity pilot, which terminated failed. |
| Protected promotion | Zero of 25 candidates; the 20 regions remain geographic reconciliations, not promoted local genetic-sharing claims. |
| Replication | No independently replicated two-trait locus. The bounded table records three exploratory exact-variant Yale labels, 14 direction-only, two no-replication, six unavailable; Yale changes the long-sleep cutoff and overlaps UKB. |
| GWAS fine-mapping | No defensible GWAS credible set or GWAS PIP for promoted regions; all 25 candidate rows are method-prerequisite holds. |
| Trait–trait coloc | No two-trait H0–H4/PP.H4 estimate. Shared coordinates or LD do not substitute. |
| Brain eQTL/sQTL coloc | No full-cis GWAS–QTL H4 estimate. Six pair candidates carry descriptive published molecular credible-set overlap only. |
| Tissue, cell type, pathways | No confirmatory enrichment estimate or causal gene assignment; source-tissue, positional and regulatory context remains descriptive. |
| Strongest supported conclusion | The atlas has validated global genetic-correlation and PLACO candidate context, with a complete bounded 25-candidate/20-region evidence package. Local shared-genetic or causal claims require the failed family gate and downstream methods to be resolved. |
| Manuscript readiness | The bounded exploratory Results/Methods/Discussion/Limitations and aligned evidence tables pass their reproducibility audit. A confirmatory manuscript is not ready. |

The final Git commit for this report is recorded in the task handoff after
the last integrity audit; source and output hashes in this report and its
linked receipts are the data identity, independently of a moving repository
HEAD.

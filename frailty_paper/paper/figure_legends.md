# Main-figure legends — provisional drafts

These drafts describe the figures currently supported by the frozen and
source-audited outputs. They are not submission-ready. Figures 3–5 remain
incomplete or unjustified as summarized below.

## Figure 1. Evidence acquisition and phenotype-specific study framework

(A) PubMed-only search snapshot: 61,009 records were identified, 4,892
duplicate occurrences were linked during deduplication, and 56,117 unique
records remain queued for screening. No records from licensed manual databases
are included, and screening has not begun; this panel is an acquisition
framework, not a completed PRISMA flow. (B) Planned staged analysis, with
primary, sensitivity, replication and downstream evidence kept distinct.
(C) The locked 12-trait sleep and circadian panel, grouped by measurement
domain. (D) Frailty constructs and current evidence status. FI and latent-factor
results are available; Fried, HFRS and physical-component sleep-pair analyses
are not currently eligible. See the source-data TSV and provenance JSON
accompanying `analysis/figure1_evidence_framework` for counts, source files and
hashes. This figure remains provisional pending licensed database exports,
screening and completion of the analyses it depicts.

## Figure 2. Genome-wide genetic correlation across sleep traits and frailty definitions

(A) LDSC genetic-correlation estimates (r_g) for the 12 locked sleep and
circadian phenotypes paired with the primary Frailty Index (FI) and seven latent
frailty factors. Cell colors encode r_g on a diverging scale clamped at −0.7 and
+0.7; cell labels are rounded to two decimal places. An asterisk marks q≤0.05
within the corresponding column's correction family: FI q-values are inherited
from the frozen atlas's all-396-pair Benjamini–Hochberg family, whereas
latent-factor q-values use the fixed 84-pair secondary family. The right panel
reports the FI family counts, latent-family count, absence of an eligible
independent pairwise sleep–frailty rg replication, and unresolved participant
overlap. Nine of 12 FI comparisons meet the inherited BH criterion and
same-family Bonferroni sensitivity; 48 of 84 latent-factor comparisons meet the
secondary-family BH criterion and 37 also meet the same-family Bonferroni
sensitivity. These families are not directly compared, and
the latent-factor results are sensitivity-only. Neither the color pattern nor
q-significance differences test which frailty dimension has a stronger
relationship with sleep. The insomnia–general-factor cell has direct
construct-level part-whole dependence because insomnia is one of the 30 input
deficits used to construct the general factor; it is not independent
corroboration. Daytime sleepiness is related but not identical to the model's
tiredness/lethargy indicator. LDSC genetic correlation does not establish
causation, local sharing or mechanism. Source estimates, standard errors, P/q
values, family labels, overlap status and provenance are in the accompanying
`analysis/figure2_global_rg_matrix_source_data.tsv` and provenance JSON.
Pair-level indicator overlap is recorded in Supplementary Table 7. This figure
is provisional and does not represent the unavailable physical, Fried
or HFRS frailty comparisons.

## Figure 3. Latent frailty-factor estimates — partial, sensitivity-only

Forest plots show LDSC genetic-correlation estimates (points) and Wald 95% CIs
(bars) for the fixed secondary family of 12 sleep/circadian traits paired with
seven latent frailty factors. Asterisks mark BH q≤0.05 within the 84-pair
family. The insomnia–general-factor pairing has direct part–whole dependence:
insomnia is one of the 30 deficits used to construct the general factor. Exact
participant overlap is unknown, and the plot does not test differences between
correlated factors or establish independent replication or causation. Physical
components, Fried phenotype and HFRS are not shown because eligible sleep-pair
estimates are unavailable. Figure 3 therefore remains partial and cannot
support a complete cross-definition frailty decomposition. Source estimates,
Wald interval derivation, input/output checksums and limitations are recorded in
`analysis/figure3_latent_factor_forest_source_data.tsv` and
`analysis/figure3_latent_factor_forest.provenance.json`.

## Figure 4. Local sharing — interim evidence-status panel

This panel reports analysis status only. The 2026-09-26 23:07:06–23:07:58 UTC
full-family receipt audit validated 18,234/29,940 receipts; inventory reached 18,237 and runner state reported 18,237 during scanning. Both workers were live on
distinct claims, with zero receipt, claim or duplicate-identity issues. All 12
trait-level failure lower bounds exceed the frozen 1% allowance, so local
estimates remain inadmissible. A failed gate does not establish absence of
local sharing. Machine-readable evidence is
`analysis/lava_full_receipt_integrity_2026-09-26_2307_checkpoint.json`,
`analysis/lava_full_receipt_failure_gate_audit_2026-09-26_2307_checkpoint.md`,
`analysis/figures4_5_downstream_evidence_status_source_data.tsv` and
`analysis/figures4_5_downstream_evidence_status.provenance.json`.

## Figure 5. Molecular follow-up — eligibility-status panel

This panel reports evidence availability, not biological findings. No eligible
frailty-specific shared-locus set is available for fine-mapping, trait–trait
colocalization or locus-matched QTL integration. eQTL Catalogue and BrainSCOPE
are resource indexes, pQTL objects remain unacquired, and muscle single-cell
resource metadata retain unresolved archive/donor/build/reuse details. No gene,
protein, tissue or cell-type claim is made. Reassess downstream molecular work
only after eligible upstream loci and locus-matched molecular data are
available. Sources, limitations and output hashes are in
`analysis/figures4_5_downstream_evidence_status_source_data.tsv` and
`analysis/figures4_5_downstream_evidence_status.provenance.json`.

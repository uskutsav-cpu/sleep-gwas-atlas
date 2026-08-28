# Novelty-Enriched Phenome Discovery Extension

This namespace is scientifically and operationally separate from the locked
`atlas-v1.0` core. The core contains 45 traits and the prespecified 396
sleep-by-non-sleep tests. Nothing under this directory changes the core panel,
its GWAS identities, its QC thresholds, or its multiple-testing family.

The extension is exploratory and will be described as:

> The prespecified 45-trait core atlas was complemented by a separately
> defined novelty-enriched phenome-wide discovery extension.

## Immutable core boundary

`core_checkpoint.json` records the producing commit plus byte-level SHA-256
hashes for the locked panel, source/schema manifests, Phase-1 QC tables, the
complete 396-pair result family, and its figure. Run the checkpoint verifier
before and after every extension stage:

```bash
python3 discovery_extension/scripts/00_verify_core_checkpoint.py
```

The verifier fails closed if any pinned core artifact changes, if the current
commit no longer descends from the checkpoint commit, or if the core result is
not the exact 12 x 33 Cartesian product.

Core analysis artifacts use exact SHA-256 verification. The shared
`environment/tool_versions.tsv` registry uses a narrower append-only policy:
every row present at the checkpoint must remain identical and in order, while
uniquely named tools for later analyses may be appended. This preserves the
recorded core runtime snapshot without treating unrelated new tooling as a
change to the locked 45-trait/396-pair analysis.

## Namespace contract

Extension-only configuration, manifests, logs, results, figures, provenance,
and acceptance gates live below this directory. Synthetic fixtures, if ever
introduced, must live under `synthetic/` and must never be copied into real
results.

## Frozen extension panel

The pre-analysis workflow is intentionally sequential:

1. `01_build_panukbb_universe.py` creates a broad 452-trait EUR universe from
   the official Pan-UKB phenotype and heritability manifests.
2. `02_build_novelty_prescreen.py` compares it with source supplements from
   four prior sleep-genetics screens. An exact no-match is explicitly not
   treated as proof of novelty.
3. `03_build_candidate_pool.py` applies the prospective source h2/intercept,
   interpretability, core-overlap, and redundancy filters, yielding 242 traits.
4. `04_lock_extension_panel.py` freezes 100 domain-stratified traits before any
   extension correlation, together with a separate 1,200-test FDR contract.
5. `05_validate_extension_panel.py` fails closed if the core checkpoint or any
   locked extension input changes.
6. `06_acquisition_preflight.py` checks whether all dense sources can be
   retained before starting a download.
7. `07_build_source_contracts.py` binds every locked trait to its exact bgzip
   and tabix object, records range-verified continuous/binary schemas, pins the
   shared variant-INFO/rsID map, and writes the prospective harmonization and
   rerun-h2 rules.
8. `08_validate_source_contracts.py` verifies those contracts against the
   panel lock and refuses to treat the source-level Pan-UKB h2 screen as the
   extension's still-pending rerun.
9. `09_build_panukbb_hm3_reference.py` performs the one-time, fail-closed
   GRCh37 coordinate/allele/rsID join between the full Pan-UKB variant manifest
   and the repository's checksum-pinned HapMap3 map.
10. `10_harmonize_panukbb.py` streams a locked source into canonical LDSC
    columns, applying the prespecified low-confidence, variant-identity,
    autosome, SNP, strand-ambiguity, INFO, MAF, MHC, and duplicate rules.

Both continuous and binary branches of the harmonizer have an isolated smoke
test at `synthetic/test_panukbb_harmonization.py`. It writes only to a temporary
directory and never creates or populates a real result path.

The downstream contract is already executable once the acquisition gate is
cleared: `11_munge_extension.sh` creates extension-only HapMap3 inputs,
`12_h2_extension.sh` reruns observed-scale h2 and applies the fixed Z/intercept
gate, `14_rg_extension.sh` runs the 12-by-pass-trait family, and
`13_collate_extension_ldsc.py` writes both the primary tested family and the
complete 1,200-pair status universe with extension-only BH FDR. The discovery
matrix records the cross-trait intercept, SE, post-merge and valid-allele SNP
overlaps, ancestry, and analysis status for each tested pair. The h2 step emits
the requested `extension_trait_readiness.tsv` plus a separate QC-failed table;
the rg step emits `extension_rg_matrix.tsv`.

`19_build_extension_views.py` creates ranked positive, ranked negative, and
domain-specific tables. `15_plot_extension.py` creates the overview heatmap,
an FDR network, and a multi-page domain-view PDF; its caption explicitly keeps
pre-screen novelty categories distinct from pair-level novelty claims. The
entire 100-trait h2, variable-size rg, FDR, SNP-diagnostic, exclusion-ledger,
ranked-view, and figure path is exercised under
`synthetic/test_extension_ldsc_collation.py` only.

`16_extension_acceptance.py` emits a machine-readable 17-stage gate table, an
adversarial-risk checklist, and a plain-language execution-status report. It
records checkpoint drift and missing downstream artifacts as blocking or
pending states instead of aborting before the failure is documented.

For result-dependent follow-up, `17_prepare_pair_novelty_audit.py` creates a
PENDING-only row for every extension-FDR-significant pair and
`18_validate_pair_novelty_audit.py` enforces the brief's six novelty classes,
requires several targeted searches for `APPARENTLY_NOVEL`, and refuses a
STRONG novelty label without successful independent replication. The canonical
artifact is `results/novelty/extension_novelty_audit.tsv`.

`20_prioritize_extension_pairs.py` applies the frozen Tier A/B/C rules across
the complete tested family and writes `novel_hit_priority.tsv`. Tier A requires
extension FDR, effect size, both h2 gates, clean analysis status, at least
500,000 valid-overlap SNPs, and Strong/Moderate `APPARENTLY_NOVEL` or
`NO_DIRECT_RG_FOUND` pair evidence. It also preserves biological plausibility,
obviousness, replication/dense-data/molecular-QTL availability, and the literal
h2/intercept values rather than hiding them in a composite score. Tier B
requires Tier A plus replication class `REPLICATED` or strong local support;
everything else remains explicitly Tier C. The replication and local/pleiotropic/fine-
mapping/colocalization/mechanistic rules are frozen in
`config/replication_contract.json` and `config/followup_contract.json`.

Replication is a second locked workflow, not a post hoc lookup:
`21_prepare_replication_queue.py` emits result-free source-curation rows for
Tier A/B pairs not yet replicated and independently locks that candidate
family; `22_lock_replication_manifest.py` requires verified local
full-resolution source SHA-256 values, completed source searches, phenotype
compatibility, EUR ancestry, and confirmed non-overlap before result access;
and `23_collate_replication.py`
classifies the exact locked family as `REPLICATED`,
`DIRECTIONALLY_CONCORDANT`, `UNDERPOWERED`, or `FAILED_REPLICATION` using the
locked Bonferroni threshold, while retaining effect heterogeneity. An
evidence-backed `NO_INDEPENDENT_DATASET` row is preserved as the explicit fifth
outcome rather than being dropped from the candidate family. Both a successful
source and an unavailable-source branch are exercised with synthetic data.

Local architecture is likewise a locked, result-preserving workflow.
`24_local_architecture_preflight.py` verifies the pinned official LAVA 0.1.5
and HDL 1.4.3 code/entrypoints independently of their LD references and dense
inputs. The source and reference manifests retain the exact Git commits,
citations, official LAVA UKB v1.1 download components, and the unresolved
checksum status rather than treating a URL as a verified file.
`25_prepare_local_analysis_queue.py` keeps Tier A/B discoveries alongside a
separate result-free curation set of globally-null pairs, so global rg
significance is not the sole entry criterion. `26_lock_local_analysis_manifest.py`
then freezes pair membership, dense-input/locus checksums, methods, and the
pair-by-locus family before local result access. Finally,
`27_collate_local_architecture.py` preserves every method/locus row, applies BH
only to the clean estimable LAVA family, keeps HDL-L as sensitivity evidence,
and emits the four required local-architecture flags plus a pair summary. The
lock and flag semantics are exercised without real data in
`synthetic/test_local_architecture_contract.py`.

The pleiotropy layer is also executable but input-gated. The official
stand-alone PLACO+ 0.2.0 source is pinned by Git commit and file SHA-256;
`28_fetch_placo_plus.sh` reproduces the checksum-verified local installation,
and `28_pleiotropy_preflight.py` verifies `var.placo`, `cor.pearson`, and
`placo.plus`. `29_prepare_pleiotropy_queue.py` freezes only independently
replicated Tier B pairs, `30_lock_pleiotropy_manifest.py` traverses and locks
the full genome-wide harmonized input plus the LD reference before result
access, and `31_run_placo_plus.R` estimates the nuisance variance/correlation
from genome-wide data and retains per-variant numerical failures.
`31_collate_pleiotropy.py` streams the complete locked result family, validates
post-scan LD clumping, and emits `novel_shared_loci.tsv` with an explicit
statistical-pleiotropy-only claim guard. The real source functions and the
end-to-end lock/runner/collator path are exercised under
`synthetic/test_pleiotropy_contract.py`.

The current acquisition preflight is deliberately fail-closed. The 100 exact
phenotype files, their indexes, and the shared Pan-UKB variant reference total
about 214.7 GiB compressed (about 246.9 GiB with the locked 1.15 safety
factor). `provenance/acquisition_preflight.json` records the volume-specific
free-space measurement and shortfall. No bulk download, harmonized statistic,
h2 rerun, genetic correlation, FDR result, or biological follow-up is claimed
until this gate passes.

The external imaging queue is deliberately not part of the 100-trait lock: its
papers and repositories are verified, but exact file identities, checksums,
schemas, and source h2/intercepts remain unresolved. Adding one later would be
a new panel version, never a replacement inside this locked analysis.

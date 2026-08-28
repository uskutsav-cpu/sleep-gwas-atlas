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
overview is produced by `15_plot_extension.py`; its caption explicitly keeps
pre-screen novelty categories distinct from pair-level novelty claims. The
entire 100-trait h2, variable-size rg, FDR, exclusion-ledger, and figure path is
exercised under `synthetic/test_extension_ldsc_collation.py` only.

`16_extension_acceptance.py` emits a machine-readable 17-stage gate table, an
adversarial-risk checklist, and a plain-language execution-status report. It
records checkpoint drift and missing downstream artifacts as blocking or
pending states instead of aborting before the failure is documented.

For result-dependent follow-up, `17_prepare_pair_novelty_audit.py` creates a
PENDING-only row for every extension-FDR-significant pair and
`18_validate_pair_novelty_audit.py` refuses incomplete reviews or a STRONG
novelty label without successful independent replication. The replication and
local/pleiotropic/fine-mapping/colocalization/mechanistic rules are frozen in
`config/replication_contract.json` and `config/followup_contract.json`.

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

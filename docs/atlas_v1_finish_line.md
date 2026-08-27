# atlas-v1.0 finish-line ledger

This ledger translates the agreed project finish line into auditable gates. It
separates completed reproducibility infrastructure from scientific results.
Pipeline code, placeholder files, and synthetic smoke-test output do not satisfy
a scientific gate.

Run the machine-readable audit with:

```bash
python3 scripts/99_atlas_acceptance.py --report-only
```

Omit `--report-only` in a release job so any incomplete gate returns a failure.

## Current state (2026-08-26)

| Gate | State | Current evidence or blocker |
|---|---|---|
| Locked scope | Complete | Exactly 12 sleep/circadian plus 33 non-sleep traits are locked by ordered-ID hash in `analysis_panel.tsv` and its lock file. |
| Reproducible Phase 0/1 architecture | Complete | Granular readiness, pinned runtimes, Snakemake entry points, CI, negative contract tests, and a full synthetic 396-pair smoke test are implemented. |
| Source curation | Blocked | 31/45 selected sources are registry-verified; 14 require a defensible selected release and file-level source record. |
| Source schemas | Blocked | `gwas_schemas.tsv` makes literal, trait-specific allele/effect/statistic mapping an explicit gate; 25/45 mappings are verified and 20 remain unresolved. IBD, Crohn, UC, and stroke are understood but explicitly require variant-to-rsID mapping before harmonization. |
| Harmonization | Blocked | Snoring has been locally materialized and harmonized (7,168,629/11,010,158 rows retained), but the ignored data artifact is absent from a fresh clone and the other 44 traits remain unresolved. |
| LDSC h2 and Phase 1 | Blocked | No real munged files or LDSC results are versioned. The synthetic outputs are tests only. |
| Full covariance and downstream science | Blocked | MiXeR, LAVA, shared-locus discovery, Genomic SEM, factor GWAS, fine-mapping, colocalization, molecular/regulatory/cell/pathway work, and MR await real Phase 0/1 inputs. |
| Integrated atlas and robustness | Blocked | Canonical evidence tables and robustness outputs do not yet exist. |
| Immutable release | Blocked | `releases/atlas-v1.0` must be produced only after every scientific gate passes. |

The immediate critical path is source curation → raw acquisition →
harmonization. The repository's larger remote experimental branch contains a
broader 149-trait scope, so it must not be merged wholesale. Evidence may be
harvested from it only when it matches one of the locked 45 selected phenotype
definitions and survives primary-source verification.

## Required gates

1. Verify all 45 selected GWAS sources, versions, phenotype definitions,
   ancestry, build, schema, effect convention, sample metadata, prevalence where
   needed, access path, and checksum.
2. Produce 45 canonical harmonized files with filter-by-filter QC ledgers and 45
   HapMap3 LDSC inputs.
3. Complete locked-criteria LDSC h2 QC for all 45, retaining explicit warnings
   and exclusions.
4. Produce the real 396-pair sleep-by-disease result family and the complete
   45-trait covariance structure required by Genomic SEM.
5. Run real univariate MiXeR before eligible bivariate MiXeR, then local
   univariate and bivariate LAVA with correction over the actual local test
   family.
6. Build one `shared_loci.tsv` from complementary PLACO and conjunction-FDR
   evidence with consistent locus definitions and effect directions.
7. Use exploratory and validated confirmatory Genomic SEM models before factor
   GWAS and Q_SNP.
8. Fine-map priority shared loci, perform multi-signal trait and molecular-QTL
   colocalization, and integrate TWAS/sQTL/pQTL/PWAS evidence.
9. Add regulatory-element, cell-type, pathway, network, and cautious
   bidirectional causal-inference layers.
10. Assemble the canonical atlas tables under `results/atlas`, complete the
    predefined robustness pass, and freeze an immutable `atlas-v1.0` release
    with checksums, configs, QC, provenance, tool versions, and commit hash.

The acceptance checker names the canonical artifact expected for each gate.
Those names are the interface between analysis modules and the integrated atlas;
changing one requires an explicit reviewed change to the checker and this
ledger.

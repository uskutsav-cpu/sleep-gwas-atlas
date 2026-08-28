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

## Current state (2026-08-27)

| Gate | State | Current evidence or blocker |
|---|---|---|
| Locked scope | Complete | Exactly 12 sleep/circadian plus 33 non-sleep traits are locked by ordered-ID hash in `analysis_panel.tsv` and its lock file. |
| Reproducible Phase 0/1 architecture | Complete | Granular readiness, pinned runtimes, Snakemake entry points, CI, negative contract tests, a checksum-pinned GRCh38-to-GRCh37 point-liftover path, and a full synthetic 396-pair smoke test are implemented. |
| Source curation | Complete | 45/45 selected sources are registry-verified. Five official FinnGen R9 Finnish/European endpoints and the Burren 2024 NFE telomere stratum explicitly replace six access/ancestry blockers while preserving locked trait IDs and recording phenotype/power differences. |
| Source schemas | Complete | All 45 literal, trait-specific allele/effect/statistic mappings are verified. The six replacement archives passed exact bytes/checksum, gzip CRC, literal-header, and complete audits totaling 115,865,922 rows; all six have pinned hg38-to-hg19 plans. |
| Harmonization | Complete | All 45 traits have real filter-by-filter harmonization ledgers and canonical HapMap3 LDSC inputs. The six large GRCh38 replacements use checksum-bound streaming HapMap3 prefilters before the registered point liftover, while their complete raw sources remain retained for locus-level work. |
| LDSC h2 and Phase 1 | Complete with exclusions | All 45 h² runs are present. Forty-three traits pass the predefined Z/intercept gate; T2D is dropped for intercept 1.315 (>1.20) and melanoma for h² Z 3.57 (<4). The exact 396 sleep×non-sleep pairs are present as 372 primary plus 24 clearly labelled QC-failed sensitivity rows. The locked-family correction yields 153 primary FDR<0.05 pairs. MS and melanoma liability h² exceed one and remain explicit interpretation/robustness warnings rather than being hidden. |
| Full covariance | Complete with warnings | The real 45-by-45 S/Rg/I matrices, 1,035-by-1,035 V matrix, and all 1,035 lower-triangle estimates pass structural validation. S has three negative eigenvalues and V is ill-conditioned; the required 1,082 jackknife blocks and HDL implementation sensitivity remain explicit downstream warnings. |
| Genomic SEM validation | Failed honestly | Odd-chromosome discovery and even-chromosome confirmation were completed for 42 QC-passing traits. None of 10 candidate models passed held-out fit plus residual-admissibility criteria; factor GWAS and Q_SNP are therefore not run. |
| Other downstream science | Blocked | Pinned LAVA, MiXeR, and PLACO+/conjFDR production contracts now fail closed on their exact missing inputs and references; real runs plus fine-mapping, colocalization, molecular/regulatory/cell/pathway work, and MR remain. |
| Integrated atlas and robustness | Blocked | Canonical evidence tables and robustness outputs do not yet exist. |
| Immutable release | Blocked | `releases/atlas-v1.0` must be produced only after every scientific gate passes. |

The immediate critical path is univariate MiXeR and local-univariate LAVA,
followed by eligible bivariate runs and the remaining prespecified downstream
modules. The complete covariance structure required by Genomic SEM is now
available and structurally validated. The repository's larger remote experimental branch contains a
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

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

## Current state (2026-08-28)

The executable acceptance audit currently passes **10/23** gates. All 107
contract tests pass, every Python file compiles, every shell script parses, and
all eight R scripts parse. The repository contains no tracked working-tree
changes at this checkpoint.

| Gate | State | Current evidence or blocker |
|---|---|---|
| Locked scope | Complete | Exactly 12 sleep/circadian plus 33 non-sleep traits are locked by ordered-ID hash in `analysis_panel.tsv` and its lock file. |
| Reproducible Phase 0/1 architecture | Complete | Granular readiness, pinned runtimes, Snakemake entry points, CI, negative contract tests, a checksum-pinned GRCh38-to-GRCh37 point-liftover path, and a full synthetic 396-pair smoke test are implemented. |
| Source curation | Complete | 45/45 selected sources are registry-verified. Five official FinnGen R9 Finnish/European endpoints and the Burren 2024 NFE telomere stratum explicitly replace six access/ancestry blockers while preserving locked trait IDs and recording phenotype/power differences. |
| Source schemas | Complete | All 45 literal, trait-specific allele/effect/statistic mappings are verified. The six replacement archives passed exact bytes/checksum, gzip CRC, literal-header, and complete audits totaling 115,865,922 rows; all six have pinned hg38-to-hg19 plans. |
| Harmonization | Complete | All 45 traits have real filter-by-filter harmonization ledgers and canonical HapMap3 LDSC inputs. The six large GRCh38 replacements use checksum-bound streaming HapMap3 prefilters before the registered point liftover, while their complete raw sources remain retained for locus-level work. |
| LDSC h2 and Phase 1 | Complete with exclusions | All 45 h² runs are present. Forty-three traits pass the predefined Z/intercept gate; T2D is dropped for intercept 1.315 (>1.20) and melanoma for h² Z 3.57 (<4). The exact 396 sleep×non-sleep pairs are present as 372 primary plus 24 clearly labelled QC-failed sensitivity rows. The locked-family correction yields 153 primary FDR<0.05 pairs. MS and melanoma liability h² exceed one and remain explicit interpretation/robustness warnings rather than being hidden. |
| Full covariance | Complete with warnings | The real 45-by-45 S/Rg/I matrices, 1,035-by-1,035 V matrix, and all 1,035 lower-triangle estimates pass structural validation. S has three negative eigenvalues and V is ill-conditioned; the required 1,082 jackknife blocks and HDL implementation sensitivity remain explicit downstream warnings. |
| Genomic SEM validation | Passed as a validated null | Odd-chromosome discovery and even-chromosome confirmation were completed for 42 QC-passing traits. None of 10 candidate models passed held-out fit plus residual-admissibility criteria; factor GWAS and Q_SNP are therefore not applicable, with immutable header-only terminal outputs rather than fabricated associations. |
| Other downstream science | Blocked | Pinned LAVA, MiXeR, PLACO+/conjFDR, fine-mapping, molecular-QTL/TWAS, regulatory, five-strategy cell-type, four-resource pathway, bidirectional-MR, atlas-edge, and robustness workflows now fail closed on exact missing inputs, releases, runtimes, or references. Only 29/45 inputs are full-resolution; the other 16 are `parkinson`, `mdd`, `ibd`, `crohn`, `uc`, `ms`, `asthma`, `t2d`, `ldl`, `hdl`, `triglycerides`, `cad`, `stroke`, `longevity`, `telomere_length`, and `melanoma`. Dense source files for 13 of those 16 remain local; the three GLGC lipid archives are absent. CATlas adult scATAC and the 16-tissue GTEx LDSC-SEG subset are source-ready (19/20 interpretation sources overall). The exact LDSC-SEG reference/runtime/trait contract is implemented, but its 1,876,474,664-byte streamed transfer has not started because current free space is below the locked 4 GiB floor. Real production runs remain. |
| Integrated atlas and robustness | Blocked | Canonical evidence tables and robustness outputs do not yet exist. |
| Immutable release | Blocked | `releases/atlas-v1.0` must be produced only after every scientific gate passes. |

## Exact production-capacity blocker

The present Mac is arm64 with 8 GiB RAM and approximately 2.2 GiB free. The
remaining primary jobs cannot be run safely on it. The fail-closed production
contracts require:

| Work family | Exact missing fixed payload | Minimum execution capacity |
|---|---:|---|
| LAVA UK Biobank LD v1.1 | 14,110,596,095 bytes compressed; 15 GiB reported extracted | 35 GiB free for download, extraction, and margin |
| MiXeR | 6,579,093,199-byte 64-file reference plus 2,106,984,629 compressed container bytes | x86_64/amd64, 32 GiB RAM, 16 physical cores recommended, 30 GiB free |
| PLACO+/conjunction FDR | 2,383,912,974-byte LD reference plus 274,423,819-byte variant template | MATLAB, SciPy, 16 GiB RAM, 20 GiB free |
| TWAS | 3,135,665,776-byte, 98-file PredictDB model family plus the pinned MetaXcan runtime | 8 GiB RAM and 20 GiB free; all 45 full-resolution GWAS |
| LDSC-SEG | 1,876,474,664 streamed bytes, of which 1,593,837,597 are temporary | 4 GiB free and explicit large-transfer/temporary-cleanup acknowledgement |
| Genome-wide GRCh37 identity map | 2,701,503,051-byte Pan-UKBB variant manifest | 10 GiB free for guarded download; 15 GiB free for the deterministic disk-backed build |
| Three missing GLGC dense GWAS | 6,844,892,917 bytes across HDL, LDL, and triglycerides | enough additional working space to harmonize and retain all 16 dense replacements |

The byte-pinned entries above total **40,013,547,124 bytes (37.266 GiB)** before
MetaXcan package payloads, dynamic QTL queries, harmonized dense outputs,
intermediate matrices, and final result artifacts. Because LAVA archives and extracted reference files coexist,
the 35 GiB single-job floor is not a safe whole-project allocation; an x86_64
host with at least 32 GiB RAM, MATLAB, and at least 100 GiB free is the
conservative production handoff target.

The immediate critical path is to build the sealed genome-wide identity map and
create the 16 full-resolution canonical GWAS on that host, then run
local-univariate LAVA and univariate MiXeR, followed by
eligible bivariate runs and the remaining prespecified modules. Fine-mapping
must wait for the complete 396-pair PLACO+/conjunction-FDR family because a
missing `shared_loci.tsv` is unknown, not a valid zero-locus result. The
complete covariance structure required by Genomic SEM is available and
structurally validated. The repository's larger remote experimental branch contains a
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
   GWAS and Q_SNP; if none validates, publish the locked not-applicable family.
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

# Track B independent-replication source audit

Frozen: 2026-09-01, before reading any new replication result.

## Pair A — snoring / parental lifespan

`NO_VALID_REPLICATION`. The exact public source is the discovery GWAS. LifeGen is an ancestor component of that meta-analysis. HRS is small, differently defined, and lacks an identified dense public release. The independent AncestryDNA component is not separately available as analysis-ready public summary statistics, while its combined meta-analysis includes UK Biobank. Extreme longevity is retained only as a nonexact sensitivity and cannot earn replication credit.

Mechanistic follow-up may proceed as an unreplicated discovery analysis, but Pair A cannot reach Track B Tier 2 until a truly independent, closely matched parental-lifespan dataset is obtained.

## Pair B — insomnia / ADHD

FinnGen R13 endpoint `F5_ADHD` is frozen as the primary external-phenotype replication source: 5,559 cases and 489,493 controls, Finnish-European, GRCh38, 802,195,177 compressed bytes. The version is pinned by GCS generation `1777989549884404` and MD5/ETag `e7aedc4602071840bbd5abe62d6fcdcc`.

This is independent by published cohort descriptions from the iPSYCH/deCODE/PGC discovery meta-analysis, not by individual-level linkage. It is much smaller in case count and uses Finnish register ascertainment, so a null result can be `UNDERPOWERED`; founder-population LD and phenotype differences must be carried into heterogeneity interpretation.

The UKB insomnia GWAS is intentionally held fixed. Track B asks for an independent external phenotype GWAS whenever feasible; replacing both GWAS would change the scientific question and reduce comparability.

## Frozen execution rule

- Use the same standardized LDSC/HapMap3 pipeline and thresholds as the atlas.
- Do not substitute a different ADHD endpoint after seeing the result.
- Record h2/QC before interpreting rg.
- A well-powered, materially incompatible opposite direction triggers the specified no-go gate.
- A low-power null is not promoted to failed replication.

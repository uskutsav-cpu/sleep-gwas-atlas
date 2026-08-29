# Novelty-Enriched Phenome Discovery Extension — final scientific report

Generated: 2026-08-29T13:12:13Z

## Publication-readiness verdict

**GLOBAL DISCOVERY AND INDEPENDENT REPLICATION COMPLETE; LOCAL AND MECHANISTIC FOLLOW-UP RESOURCE-BLOCKED.**

The immutable 45-trait/396-pair core remains separate and unchanged. This separately locked extension prospectively selected 100 phenotypes, executed all 1,200 planned sleep-by-phenotype LDSC tests, audited every extension-FDR hit, and preserved the complete 217-pair replication candidate family. Findings are reportable as replicated global genetic correlations with explicit novelty and causal-language limits; they are not local-sharing, pleiotropy, colocalization, or mechanism claims.

## Required scientific summary

| # | Required count | Value | Status |
|---:|---|---:|---|
| 1 | candidates considered | 242 | COMPLETE |
| 2 | selected before analysis | 100 | COMPLETE |
| 3 | passing rerun h2 qc | 100 | COMPLETE |
| 4 | sleep x trait tests executed | 1200 | COMPLETE |
| 5 | extension fdr significant pairs | 603 | COMPLETE |
| 6 | already established pairs | 1 | COMPLETE |
| 7 | partial extension pairs | 294 | COMPLETE |
| 8 | no direct prior rg found pairs | 308 | COMPLETE |
| 9 | apparently novel pairs | 0 | COMPLETE |
| 10 | independently replicated pairs | 23 | COMPLETE |
| 11 | pairs with significant local correlations | NA_BLOCKED_UPSTREAM | BLOCKED_UPSTREAM |
| 12 | pairs with pleiotropic loci | NA_BLOCKED_UPSTREAM | BLOCKED_UPSTREAM |
| 13 | prior robust colocalized signals | NA_BLOCKED_UPSTREAM | BLOCKED_UPSTREAM |
| 14 | pairs with candidate genes or mechanisms | NA_BLOCKED_UPSTREAM | BLOCKED_UPSTREAM |

## Discovery and novelty results

All 100/100 extension traits passed the prespecified rerun h2 gates. The exact 1200-pair family yielded 603 BH-FDR-significant pairs; 381 also had |rg|>=0.15. The maximum absolute cross-trait LDSC intercept was 0.1386, and the minimum valid-allele overlap was 1,128,043.

The literature audit covered 603/603 FDR-significant pairs: 1 already established or same-pair/new-dataset result, 294 partial extensions, 308 with no direct prior rg found, and 0 classed APPARENTLY_NOVEL. The zero APPARENTLY_NOVEL count is deliberate: none cleared the audit's strongest novelty bar, and no first-ever claim is made.

## Independent replication

Among the frozen 217 candidates, 23 met the 0.05/217 replication threshold, 18 more were concordant but below that threshold, 17 were stopped by prespecified replication h2/intercept gates, and 159 had no suitable independent dataset. None of the 41 tested pairs was direction-discordant. Of the 23 replicated pairs, 18 used exact phenotype matches and 5 used documented comparable definitions; 7 showed effect-size heterogeneity P<0.05 and are retained with that qualification.

The replicated set concentrates in gastrointestinal, pulmonary, smoking-related, pain, and musculoskeletal phenotypes. This pattern is biologically interpretable but cannot distinguish direct biology, vertical or horizontal pleiotropy, mediation, ascertainment, or residual sample-structure effects.

## Top replicated underreported discoveries

| Rank | Sleep trait | External phenotype | Discovery rg | Replication rg | Replication P | Literature status |
|---:|---|---|---:|---:|---:|---|
| 1 | insomnia | Abdominal pain | 0.559 | 0.410 | 5.78e-35 | NO_DIRECT_RG_FOUND |
| 2 | insomnia | K21 Gastro-oesophageal reflux disease | 0.522 | 0.364 | 3.11e-22 | NO_DIRECT_RG_FOUND |
| 3 | insomnia | Diaphragmatic hernia | 0.432 | 0.192 | 9.34e-07 | NO_DIRECT_RG_FOUND |
| 4 | shortsleep | Abdominal pain | 0.407 | 0.318 | 4.25e-18 | NO_DIRECT_RG_FOUND |
| 5 | insomnia | Noninfectious gastroenteritis | 0.521 | 0.306 | 3.73e-07 | NO_DIRECT_RG_FOUND |
| 6 | insomnia | Tobacco use disorder | 0.387 | 0.426 | 1.49e-11 | NO_DIRECT_RG_FOUND |
| 7 | insomnia | J44 Other chronic obstructive pulmonary disease | 0.435 | 0.326 | 5.1e-23 | NO_DIRECT_RG_FOUND |
| 8 | shortsleep | K21 Gastro-oesophageal reflux disease | 0.382 | 0.196 | 4.06e-08 | NO_DIRECT_RG_FOUND |
| 9 | shortsleep | Tobacco use disorder | 0.441 | 0.361 | 4.76e-08 | NO_DIRECT_RG_FOUND |
| 10 | shortsleep | Noninfectious gastroenteritis | 0.424 | 0.206 | 0.000111 | NO_DIRECT_RG_FOUND |

The complete 23-row ranked table is `results/top_novel_discoveries.tsv`. Local loci, pleiotropic loci, colocalized genes, and mechanistic evidence are explicitly `NA_BLOCKED_UPSTREAM`, not zero.

The complete replication scatter, effect-estimate forest plot, and locked-family attrition flow are in `figures/extension_replication_summary.png` and its vector PDF counterpart.

## Resource-bounded downstream follow-up

The original full Pan-UKB mirror would require 246.91 GiB under the locked safety factor. Streaming and receipt sealing solved the LDSC acquisition problem for all 100 traits, but did not create the dense local-analysis inputs. At the final preflight only 1.280 GiB remained. The published LAVA UKB v1.1 reference alone is 15 GiB unzipped; the signed fine-mapping LD and PLACO+ clumping reference are also absent. The retained HapMap3-munged LDSC files omit contract-required dense locus coverage and per-variant MAF/INFO and therefore were not silently relabeled as full-resolution inputs.

A result-free local queue preserves 217 priority discoveries and 597 globally-null secondary candidates. A separately locked PLACO+ candidate family preserves all 23 replicated Tier-B pairs. No LAVA, HDL-L, PLACO+, SuSiE, coloc, molecular-QTL, or mechanism result is reported.

## Claim boundary

The strongest supported statement is that 23 underreported sleep-phenotype genetic correlations replicated directionally in independent EUR GWAS at the frozen family-wise threshold. Genetic correlation is not causation. No result is described as first-ever, locally shared, pleiotropic, colocalized, fine-mapped, gene-mediated, or mechanistic. Seven replicated effects require explicit heterogeneity disclosure, and five use comparable rather than exact phenotype definitions.

# Methods map: CDG3 → this project

Reference: Grotzinger, Werme, Peyrot et al., "Mapping the genetic landscape
across 14 psychiatric disorders," *Nature* 649:406–415 (2026).
Companion: Romero et al., "Exploring the genetic overlap between twelve
psychiatric disorders," *Nat Genet* (2022) — doi 10.1038/s41588-022-01245-2.

Print this. It is the one page that converts the reading into work, and it
doubles as a meeting handout.

| Their figure | What it shows | Tool | Input it needed | Our equivalent | Status |
|---|---|---|---|---|---|
| Fig. 1a (lower) | rg heatmap, 14 disorders | LDSC bivariate | munged HapMap3 sumstats + EUR LD scores | Fig 2: 6 sleep × N disease | `04_rg.sh` → `06_heatmap.py` |
| Fig. 1a (upper) | correlations among factors | GenomicSEM | LDSC genetic covariance matrix | after QC gate decides factor count | Phase 4 |
| Fig. 1b | five-factor model diagram | GenomicSEM EFA→CFA | same covariance matrix | Sleep Disturbance / Circadian factors — **may not both survive** | Phase 4 |
| Fig. 1c | hierarchical *p*-factor | GenomicSEM 2nd-order | same | a general "sleep axis" factor? open question | Phase 4 |
| Fig. 2 | rg of factors vs 31 external traits | GenomicSEM + Q_Trait | factor model + external sumstats | sleep factors vs disease domains | Phase 4/5 |
| Fig. 3a | local rg heatmap + network | LAVA v0.1.0 | sumstats + 1000G EUR + LD blocks + LDSC intercepts | top pairs only | Phase 2 |
| Fig. 3b | chr11 hotspot detail | LAVA | same | — | Phase 2 |
| Fig. 4a | CC-GWAS loci heatmap | CC-GWAS | sumstats + prevalences | probably skip (needs case/control pairs) | — |
| Fig. 4b,c | Miami + QQ for factors | GenomicSEM `userGWAS` | SNP-level sumstats, all traits | factor GWAS | Phase 4.4 |
| Fig. 5a | GO enrichment | ToppGene / MAGMA | target gene lists | Phase 9 | later |
| Fig. 5c,d | cell-type enrichment | EWCE + MAGMA | snRNA-seq references | **our claimed innovation** (Phase 8) | later |
| — | stratified enrichment | stratified GenomicSEM | 162 functional annotations | Phase 8 | later |

## Parameters lifted directly from their Methods

Cite these when asked why. They are in `01_harmonize.py` and `03_h2_qc.sh`.

- Strand-ambiguous SNPs removed
- INFO > 0.6; MAF > 1%
- MHC region excluded from all sumstats before analysis
- SNPs kept where the SNP-specific effective N is >50% of the total
- Genome build GRCh37/hg19 throughout
- munge restricts to HapMap3 and INFO > 0.9
- Liability-threshold model for case/control traits; ascertainment correction
  uses the **sum of effective sample sizes**, not raw N
- EFA on even chromosomes → CFA on odd chromosomes, to avoid overfitting
- Factor loadings > 0.3 to assign an indicator to a factor
- Model fit judged on CFI and SRMR (their five-factor: CFI 0.971, SRMR 0.063)

## The two things worth raising with your mentor

1. **MiXeR is in their pipeline and absent from our proposal.** They used it
   because LDSC underestimates overlap when shared variants have mixed effect
   directions. Sleep traits plausibly have exactly that structure with
   metabolic disease. Worth asking whether to add it.

2. **Their MiXeR inclusion threshold was N_eff × h² > 12,000.** Apply that
   arithmetic to the sleep traits before promising a two-factor sleep model.
   Several will not clear it.

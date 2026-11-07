# MAD1L1 methylation instruments and sleep–frailty candidate: reference LD

[Qian et al.](https://www.nature.com/articles/s44220-026-00660-2) list three genetic instruments for their blood methylation–ADHD analysis of **cg19624444**: rs56305291, rs1637744 and rs1108036 (their Supplementary Table S2a, row 4). The separately selected sleep–frailty candidate **rs28528997** is at GRCh38 chr7:2,034,487. The paper's supplementary table places the CpG at chr7:2,008,862–2,008,863, so the candidate is approximately 25.6 kb away. The candidate's position and allele representation remain subject to the original sleep–frailty source's build and multi-allelic caveats.

To move beyond proximity, we checked phased LD in 503 European-labeled 1000 Genomes Phase 3 reference participants. Each row counts the G allele in the biallelic reference representation; the published methylation effect alleles have not been established here.

| Published methylation instrument | Reference r² with rs28528997 | Bootstrap 95% interval for r² |
|---|---:|---:|
| rs56305291 | **0.914** | 0.870–0.956 |
| rs1637744 | 0.016 | 0.003–0.039 |
| rs1108036 | 0.012 | 0.002–0.032 |

The candidate therefore tags one published methylation instrument strongly in this reference panel; the other two are much less correlated with it. This supplies a concrete variant-level reason to investigate the *MAD1L1* region across the projects. It **does not establish** that rs28528997 changes methylation, that the methylation association and sleep/frailty associations have the same causal variant, that *MAD1L1* mediates any phenotype, or that ADHD and sleep share a causal signal. The confidence intervals describe resampling of this small reference panel under an independence assumption; they are not intervals for study-specific LD. Both rs28528997 and rs1637744 have additional alleles in Ensembl beyond the biallelic 1000 Genomes records, so source allele normalization must be confirmed. No methylation effect allele was inferred from the sign of LD.

Qian et al.'s published MR uses multiple instruments and assumptions, and its molecular sample is adolescent blood. The selected sleep–frailty candidate was examined after earlier findings, its source build is not independently authenticated, and its MVP insomnia/FinnGen HFRS comparison is not independently validated for both original FI and sleep phenotypes. The Qian PGC ADHD GWAS is not a new independent replication of the original ADHD association. Definitive follow-up requires dense and suitably licensed methylation QTL and sleep/ADHD/frailty summary statistics, ancestry-matched LD, allele harmonization, multi-signal colocalization and matched-cell experiments.

[Protocol](PROTOCOL.md), [all pairwise results](instrument_ld.tsv), [script](reference_ld.py), [reference receipt](instrument_ld_receipt.json), [publisher supplement](../../PUBLICATION.md#external-dependencies) (external target: `qian_supplement.xlsx`) and saved Ensembl responses preserve the audit. The VCF was accessed by indexed remote queries; its full bytes were not hashed. The supplementary workbook and article are [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); cite Qian et al., Nature Mental Health (2026), DOI 10.1038/s44220-026-00660-2.

A subsequent [exact-variant GTEx lookup](gtex/RESULTS.md) returned eQTL associations for an antisense transcript to *MAD1L1* in nucleus accumbens and for *MRM2* in caudate for both linked variants. It also returned MAD1L1 splice associations outside brain. These linked-marker observations share donors and are not independent validation or disease–QTL colocalization.

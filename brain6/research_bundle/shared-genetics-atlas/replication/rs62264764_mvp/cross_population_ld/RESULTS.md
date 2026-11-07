# rs62264764 and rs62264767 across reference groups

The two chromosome 3 variants remain **strongly linked in all five examined 1000 Genomes Phase 3 reference groups**. Their counted alleles are A at rs62264764 and C at rs62264767. This makes it difficult to separate their genetic effects using these reference data alone, despite the different AlphaGenome predictions reported in the [earlier analysis](../RESULTS.md).

| Reference group | People with both phased calls | rs62264764 A frequency | rs62264767 C frequency | Haplotype r² (bootstrap 95% interval) |
| --- | ---: | ---: | ---: | ---: |
| African-labeled (AFR) | 661 | 0.128 | 0.129 | 0.987 (0.964–1.000) |
| Admixed American-labeled (AMR) | 347 | 0.079 | 0.079 | 1.000 (1.000–1.000) |
| East Asian-labeled (EAS) | 504 | 0.060 | 0.060 | 1.000 (1.000–1.000) |
| European-labeled (EUR) | 503 | 0.164 | 0.164 | 1.000 (1.000–1.000) |
| South Asian-labeled (SAS) | 489 | 0.125 | 0.126 | 0.991 (0.969–1.000) |

The intervals resample individuals within each reference group. When all observed haplotypes match, the bootstrap interval is degenerate at 1; it does **not** establish perfect linkage in a wider population. Two of 1,322 AFR haplotypes and one of 978 SAS haplotypes had the second counted allele without the first; none had the first without the second. Phase is inferred in this reference and group labels are broad, so the results do not characterize every ancestry or individual. r² estimates, frequencies, imputation and phase in the actual GWAS cohorts may differ.

This is **reference LD**, not replication of insomnia, ADHD or a shared causal variant. A different pattern in another cohort might improve fine-mapping, but these data do not demonstrate such resolution. Neither the AlphaGenome prediction nor this LD table identifies LSAMP or LINC03051 as the disease effector; a measured allele-specific perturbation and suitably powered, independently characterized disease cohorts remain needed. The existing local-QC failures are unchanged.

The [pre-result plan](PROTOCOL.md), [code](analyze.py), [complete counts and intervals](group_ld.tsv), [reference receipt](receipt.json), [independent count-based verification](verification.json), and [run log](run.log) preserve the check. Both exact GRCh37 coordinates and allele contrasts were matched in the [UCSC-hosted reference](https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/), and the European estimate reproduces the earlier result. The VCF and index were accessed remotely and not whole-file hashed; the local sample-panel hash is in the receipt. The verification checks group coverage and recomputes frequencies and r² from aggregate counts, but does not independently verify the downloaded phase or regenerate bootstrap samples. No individual genotypes were saved.

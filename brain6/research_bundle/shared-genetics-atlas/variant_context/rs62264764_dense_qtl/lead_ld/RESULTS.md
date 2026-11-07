# Chromosome 3 lead variants: reference-panel LD

The strongest adrenal *LSAMP* eQTL in the previously matched region is **rs6804394** (nominal eQTL p=3.94×10⁻³⁴). The smallest beta/SE-derived Brain6 insomnia p-value is at **rs1456193** (2.13×10⁻⁸); the corresponding ADHD lead is **rs73163723** (2.38×10⁻⁷). These are post-selected leads within the matched 1-Mb window and are not necessarily causal variants.

In 503 European-labeled 1000 Genomes Phase 3 reference participants, their pairwise linkage with rs62264764 and each other was:

| Pair | Reference r² | Bootstrap 95% interval |
|---|---:|---:|
| adrenal eQTL lead rs6804394 – rs62264764 | 0.061 | 0.037–0.091 |
| adrenal eQTL lead rs6804394 – insomnia lead rs1456193 | 0.026 | 0.010–0.049 |
| adrenal eQTL lead rs6804394 – ADHD lead rs73163723 | 0.300 | 0.263–0.340 |
| insomnia lead rs1456193 – rs62264764 | 0.821 | 0.758–0.878 |
| ADHD lead rs73163723 – rs62264764 | 0.316 | 0.248–0.385 |

The low reference LD between the strongest adrenal eQTL and the highlighted rs62264764/insomnia markers is consistent with why the earlier single-signal model favored different regional signals. It does **not** prove that there are distinct causal variants or exclude a weaker shared signal; single-marker LD is not conditional association or multi-signal colocalization. The highest-PIP or causal variant need not be the lead p-value SNP. The reference does not contain the original study outcomes, may not match their populations, and the bootstrap intervals reflect reference-panel sampling under an independence assumption. Signed r values use the listed alternate alleles and must not be interpreted as disease-risk directions.

[Protocol](PROTOCOL.md), [all six pairs](pairwise_ld.tsv), [code](analyze.py), and [receipt](receipt.json) document the post-hoc selection, exact allele mapping and reference hashes. The remote VCF was queried by index and not fully hashed. The underlying [dense QTL comparison](../RESULTS.md) remains conditional on its own prior, sample and variant-coverage assumptions; none of these checks resolves the failed LAVA gate.

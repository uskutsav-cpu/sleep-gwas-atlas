# European reference population holdouts

We repeated the post-selected FI univariate calculation after omitting each of the five 1000 Genomes European panel populations in turn. All five subsets retained all 2,738 matched markers, and each produced a positive FI variance diagonal. The estimates ranged from **+1.959×10⁻⁶ to +2.710×10⁻⁶**, compared with **+1.881×10⁻⁶** in the full 503-person panel. The subsets contained 396–412 people and retained 159–171 LD components at the 99% rule, versus 183 in the full panel.

| Population omitted | Remaining people | Retained components | FI variance diagonal |
|---|---:|---:|---:|
| CEU | 404 | 171 | +2.303×10⁻⁶ |
| FIN | 404 | 169 | +1.959×10⁻⁶ |
| GBR | 412 | 166 | +2.505×10⁻⁶ |
| IBS | 396 | 161 | +2.710×10⁻⁶ |
| TSI | 396 | 159 | +2.615×10⁻⁶ |

These labels are those in the reference panel. The findings show that no *single* omitted population was necessary for a positive estimate in this set of calculations. The subsets overlap heavily, are smaller than the full panel and retain different numbers of LD components, so the differences mix population composition with sample size and cannot be interpreted as ancestry effects or independent replications. They do not determine the correct LD for the GWAS or validate the estimated variance. No bivariate local correlation was computed and the original LAVA family gates remain failed.

[Protocol](PROTOCOL.md), [code](run.py), [all results](results.json), and [verification](verify.py) preserve the test. The code reconstructs public-reference genotypes in memory and stops before the parent script's output-writing section; no individual genotypes or sample identifiers are saved. The verifier checks group counts and numerical arithmetic, not sample independence, the remote VCF/index integrity or statistical calibration. The [parent sensitivity](../RESULTS.md) records the reference and allele-alignment limitations.

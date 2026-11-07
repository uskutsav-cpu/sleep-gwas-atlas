# rs1456193 reporter results in additional cell lines

The published Siraj 2026 OL33 reporter library contains **rs1456193** in all four additional cell-line files checked. None of the source allele-skew comparisons meets its reported FDR<0.05 threshold; neither allele meets the source Bonferroni activity threshold in these rows. These are context-specific measurements, not proof of no regulatory effect.

| Cell line | Reported log₂ allele skew (SE) | Nominal p | Source skew FDR |
| --- | ---: | ---: | ---: |
| GM12878 | +0.084 (0.270) | 0.754 | 0.955 |
| HepG2 | +0.351 (0.259) | 0.176 | 0.552 |
| A549 | −0.057 (0.138) | 0.680 | 0.930 |
| K562 | −0.300 (0.174) | 0.084 | 0.302 |
| SK-N-SH, previously checked | +0.00094 (0.129) | 0.994 | 0.999 |

The four new files represent lymphoblastoid, liver-derived, lung-derived and leukemia-derived cell lines; SK-N-SH is neural-derived but is not mature primary brain. A reporter construct measures activity outside its endogenous genomic environment. We do not infer disease-risk direction from the assay skew or call an expression-modulating variant from skew alone. The source FDR is the authors' adjustment within their assay family and does not account for the wider post-hoc Brain6 search.

rs1456193 is linked to the user-prioritized **rs62264764** at r²≈0.821 in the earlier 503-person European reference; it is a different variant. Its nonsignificant reporter results do not test rs62264764 directly, rule out an effect in another cell type or developmental stage, identify LSAMP/LINC03051 as a target, validate or refute AlphaGenome's accessibility prediction, or replicate a sleep–ADHD association. Reference LD and the association datasets have separate limitations described in the [prior candidate review](../../linc03051_laddr_2026/RESULTS.md).

The [pre-result plan](PROTOCOL.md), [query code](query.py), [source rows](source_rows.tsv), [member receipts](receipt.json), [verification](verification.json), and [run log](run.log) preserve the lookup. Each selected archive member passed ZIP CRC and gzip reading; the complete 929-MB archive was not downloaded or independently hashed. P-values and source FDRs are transformed from the supplied −log₁₀ columns, and exact position/alleles were matched using the previously established GRCh37 conversion. The [earlier SK-N-SH result](../RESULTS.md) came from the same published assay library, not an independent experiment.

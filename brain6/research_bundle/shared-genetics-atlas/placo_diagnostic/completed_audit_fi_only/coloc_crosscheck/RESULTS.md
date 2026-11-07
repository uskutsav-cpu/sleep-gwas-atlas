# Chr7: FI diagnostic markers versus FI–insomnia coloc candidates

The FI PLACO+ diagnostic lead on chr7, rs2396766, and the previously published Song candidate rs12671330 were compared with the five highest **conditional H4** variants from the separate FI–FinnGen-insomnia single-signal colocalization. All seven variants had unique matching GRCh37 reference records; the ten requested comparisons used phased data from 503 1000 Genomes EUR individuals.

| Conditional H4 marker | Conditional probability in FI–insomnia coloc | r² with either diagnostic anchor | 1,000-resample individual-bootstrap r² interval |
|---|---:|---:|---:|
| rs10953765 | 0.2154 | 0.714 | 0.662–0.763 |
| rs1456031 | 0.2151 | 0.747 | 0.696–0.794 |
| rs936145 | 0.0342 | 0.996 | 0.988–1.000 |
| rs12154391 | 0.0303 | 1.000 | approximately 1.000–1.000 |
| rs4730637 | 0.0296 | 1.000 | approximately 1.000–1.000 |

Both diagnostic anchors have the same listed r² values in this reference panel; their reference r² with each other is approximately 1. The table is consistent with the two analyses tagging related regional signals, but it is not independent evidence: both reuse the same FI and registry-insomnia GWAS, and the anchors and top coloc variants were selected after seeing results. The two leading probabilities are **conditional on H4**, not unconditional causal probabilities. The chr7 coloc H4 changes from about 0.587 to 0.125 as p12 changes from 10⁻⁵ to 10⁻⁶, and the PLACO diagnostic has unresolved calibration concerns. Neither analysis resolves multiple causal signals or establishes that insomnia and frailty share a causal variant.

The 503-person reference includes Finnish and other European participants but is not either study cohort. Bootstrap intervals represent resampling of reference individuals and do not capture source-population mismatch. This evidence does not assign FOXP2 or any other target gene, establish sleep–ADHD–frailty three-trait causality, or repair the failed local-correlation gate. The earlier [reference comparison](../chr7_ld/RESULTS.md) found only r²≈0.08–0.11 between these diagnostic markers and two other selected regional anchors, illustrating why broad physical-window overlap alone is insufficient.

`pairs.tsv` contains counted alleles, signed r, r² and reference intervals; `receipt.json` records source and panel hashes; `calculate.py` and the [protocol](PROTOCOL.md) reproduce the analysis. See the [FI diagnostic](../RESULTS.md) and [regional coloc](../../../frailty_coloc/RESULTS.md) for source and model limitations.

A separate [allelic sign check](direction/RESULTS.md) finds that aligned FI and registry-insomnia effect signs are compatible with reference LD for the four diagnostic-anchor/leading-H4 pairs. This is not independent causal evidence.

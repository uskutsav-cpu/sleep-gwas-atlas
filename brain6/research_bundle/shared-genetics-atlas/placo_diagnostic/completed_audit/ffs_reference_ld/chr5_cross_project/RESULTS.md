# Chromosome 5: connecting selected FFS, FI and Brain6 markers

Two of the 27 post-selected insomnia–FFS PLACO+ numerical crossings, rs1363103 and rs13164188, are the top two conditional variants in the earlier FI–registry-insomnia coloc model. Their identity is an observation about selected markers, **not** independent evidence that all traits share a causal variant: the analyses reuse insomnia data, the frailty samples overlap, and the FI–insomnia H4 posterior is prior-sensitive (0.512 at p12=10⁻⁵; 0.095 at p12=10⁻⁶). Rs1363103 was previously reported as an FI lead. In the FFS diagnostic their respective source insomnia p-values are 7.34×10⁻⁴ and 7.47×10⁻⁴, and FFS p-values 9.1×10⁻¹³ and 8.0×10⁻¹³; neither is a corrected max-p conjunction discovery. Both PLACO values remain below the numerical screening cutoff in the p-derived and beta/SE sensitivities.

The table compares selected variants in the 1000 Genomes Phase 3 reference. r² describes association between allele indicators in the reference panel, not correlation between GWAS effects. The Finnish subset is nested within the EUR panel and has only 99 individuals.

| Post-selected FFS marker | Earlier anchor | EUR r² (503 people) | FIN r² (99 people) | Interpretation |
|---|---|---:|---:|---|
| rs6596578 | rs7703746, published FFS lead | 0.957 | 0.980 | Strong tagging in these reference samples; FFS association is not a new independent FFS discovery. |
| rs6596578 | rs1363103, FI–insomnia conditional leader | 0.549 | 0.494 | Moderate linkage. |
| rs6596578 | rs2431108, Brain6 insomnia–ADHD candidate | 0.519 | 0.405 | Moderate linkage; does not imply a shared three-trait causal signal. |
| rs60271 | rs2431108, Brain6 candidate | 0.901 | 0.918 | Strong reference linkage; selected after inspecting this region. |
| rs1363103 | rs13164188, FI–insomnia conditional leaders | 0.988 | 1.000 | They tag very similar variation in this panel and are not independent observations. |

The rs60271 row is a hypothesis for follow-up, not a validated sleep–ADHD–FFS locus. It is already among the Brain6 AlphaGenome pilot variants: insomnia PIP 0.0414, ADHD PIP 0.00255, and local adult frontal-cortex DNase ALT-versus-REF log₂ ratio +0.0120 in the 501-bp window. The PIPs come from the Brain6 fine-mapping model and the DNase value is a prediction, not an experimentally measured effect or a gene assignment. For rs60271 the source insomnia p=1.72×10⁻⁴ and FFS p=1.2×10⁻⁹; its diagnostic PLACO+ p=2.54×10⁻⁹ is sensitive to the method's unresolved calibration, not confirmation of pleiotropy.

The PLACO screen was restricted to HapMap3 and its one-trait-only null stress tests showed elevated nominal rejection; stringent-tail calibration is unverified. FI includes sleep items, FFS includes exhaustion, FFS effect orientation is supported by a later 37-lead source audit but its complete source adapter remains unverified, and both insomnia sources and frailty cohorts have overlapping participants. FIN estimates are noisy and not an independent replication. Pairwise LD and separate two-trait analyses cannot resolve distinct signals, assign a gene, or establish a common causal variant across three or more traits. Appropriate conditional fine-mapping, matched LD, component-excluded phenotypes and independent cohorts remain necessary.

The [protocol](PROTOCOL.md), `pairs.tsv` and `receipt.json` preserve all 216 comparisons, counted alleles, missing reference matches, sample counts and hashes. The wider [FFS LD analysis](../RESULTS.md), [FI–insomnia coloc](../../../../frailty_coloc_chr5/RESULTS.md), and [FFS prior study review](../../../../replication/ye_2023/REVIEW.md) contain the separate evidence.

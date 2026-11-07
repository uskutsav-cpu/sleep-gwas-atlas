# Reference LD for selected FFS diagnostic markers

We compared the 66 post-selected insomnia–FFS PLACO+ threshold-crossing markers in the 1000 Genomes Phase 3 broad EUR phased reference (503 people). The table gives r² to the lowest-PLACO-p selected marker on each chromosome, among the other markers with matching reference records. The graph links selected markers with pairwise r² ≥ 0.8; component sizes are shown to describe clustering, **not** to count independent signals. A connected graph can include pairs that are not strongly linked to each other.

| Chr | Selected | Exact reference matches | Lowest-p anchor | r² to anchor range | r² ≥ 0.8 graph component sizes |
|---|---:|---:|---|---:|---|
| 5 | 27 | 24 | rs6596578 | 0.468–0.949 | 21, 2, 1 |
| 11 | 9 | 9 | rs2847308 | 0.923–0.995 | 9 |
| 16 | 4 | 4 | rs6565300 | 0.985–1.000 | 4 |
| 17 | 14 | 14 | rs12602912 | 0.988–1.000 | 14 |
| 22 | 12 | 12 | rs9611602 | 0.366–1.000 | 11, 1 |

Three chromosome 5 markers (rs4703307, rs10073929 and rs40465) lacked an exact rsID, position and allele match in this reference and are not placed in the LD graph. Their absence does not prove they are invalid or independent. Allele pair, eligible sample counts and all pairwise values are retained in the receipts and TSVs.

On chr22 the FFS marker rs5995992 has r²=0.890 with FI marker rs4820425 and 0.886 with FI marker rs926914 (502 and 503 reference individuals respectively). However, rs5995992's FFS diagnostic p moves above the numerical threshold in the beta/SE sensitivity. The lowest-p FFS anchor rs9611602 has only r²=0.317–0.323 with those two FI markers. These post-selected comparisons show potentially related regional signals and also heterogeneity within the interval; they do not establish one causal variant or replication across traits. FI and FFS samples overlap.

The broad EUR panel is not study-specific Finnish or UKB/TwinGene LD. Missing, unphased, and third-allele calls are excluded pairwise. Proximity selection and graph threshold are post hoc. The PLACO null stress and stringent-tail calibration remain unresolved; neither the number of graph components nor high r² validates discoveries, assigns genes, replaces conditional fine-mapping, or establishes sleep–frailty causality. See the [full diagnostic](../RESULTS.md) and [protocol](PROTOCOL.md).

A [chromosome 5 cross-project follow-up](chr5_cross_project/RESULTS.md) compares FFS markers with prior FI, FFS and Brain6 anchors in EUR and FIN reference samples. Two selected FFS markers are identical to the FI–insomnia conditional leaders; these analyses reuse data and do not independently establish shared causality.

# Brain-cell NUDT12 expression lookup for chromosome 5 reporter variants

**The available astrocyte data do not support a detectable NUDT12 expression association for either selected reporter variant.** We queried the complete NG00184 nominal single-nucleus eQTL chromosome 5 files for rs6874138 and rs325485 in astrocytes and oligodendrocyte progenitors. Two astrocyte records matched the exact GRCh38 position, alleles and NUDT12 target; neither matched record was returned in the oligodendrocyte-progenitor file.

| Variant | Astrocyte ALT-aligned beta (SE) | Nominal p | Source-reported target FDR | Oligodendrocyte progenitors |
|---|---:|---:|---:|---|
| rs6874138 T>C | −0.00666 (0.03646) | 0.855 | 0.980 | No matching row |
| rs325485 A>G | −0.00658 (0.03774) | 0.862 | 0.980 | No matching row |

Both nominal p-values exceed 0.0125, the Bonferroni threshold for the four selected variant–cell slots. This small post-selected family is not a correction for all earlier genomic and functional searches. Source FDR is the source's target-level statistic and should not be treated as a correction across Brain6's candidate selection. The [portal](https://xqtl.niagads.org/download) returns HMT-filtered associations; its per-rsID downloads showed no records for either variant. Missing rows and filtered portal responses do not prove no association or no test.

The astrocyte records report each variant as outside the standard 1-Mb cis window but within the source's TAD-enhanced window: rs6874138 is 1,001,106 bases from the annotated NUDT12 TSS and rs325485 is 1,096,878 bases away. The absence of a progenitor row may reflect source windowing or other filtering; its cause is not established. Both cell datasets derive from an aging-brain collection and may share donors. The effect estimates are not clinical risk effects, and the correct effective sample size and source-model assumptions were not independently verified. These observations do not rule out effects in another cell type, developmental state, molecular phenotype or tissue, nor do they validate or falsify the authors' ABC model.

The [protocol](PROTOCOL.md), [exact records](exact_qtl_rows.tsv), [retrieval and hash receipt](receipt.json), [code](lookup.py) and [verification](verification.json) document the four slots. Both full source members and indices match the previously recorded publisher MD5s; tabix warned that index timestamps were older than the files, but the recorded index and data hashes matched. Checks verify extracted rows and saved responses, not statistical calibration or biological causality. Data were prepared and distributed by the National Institute on Aging Alzheimer's Disease Data Storage Site (NIAGADS), University of Pennsylvania (U24-AG041689); accession **NG00184**.

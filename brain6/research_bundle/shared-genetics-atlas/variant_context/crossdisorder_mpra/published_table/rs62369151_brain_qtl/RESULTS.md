# Brain-cell expression context for rs62369151

The reporter candidate **rs62369151** and its linked ADHD-prioritized neighbor **rs4916723** both have measured expression associations with the two previously selected brain-cell targets. For each variant, the source reports an association with LINC02060 in astrocytes and TMEM161B-AS1 in oligodendrocyte progenitors. Both variants' ALT alleles are C; the following are source ALT-aligned betas, not fold changes or disease-risk effects.

| Variant | Astrocyte LINC02060 beta (SE), p | Progenitor TMEM161B-AS1 beta (SE), p |
|---|---|---|
| rs62369151 | 0.1669 (0.02531), 8.78×10⁻¹¹ | −0.1709 (0.02726), 6.50×10⁻¹⁰ |
| rs4916723 | 0.1812 (0.02473), 6.99×10⁻¹³ | −0.1779 (0.02655), 4.45×10⁻¹¹ |

All four nominal p-values are below 0.0125, the Bonferroni threshold for these four **post-selected** variant–cell slots; this does not correct the full set of earlier candidate, gene and cell searches. Source-reported target FDRs are 1.76×10⁻⁹ and 2.02×10⁻⁸ for the two rs62369151 rows, which also appear in the HMT-filtered portal. The two markers are linked (phased r²≈0.833 in 503 European-labeled reference individuals), and the cell analyses may share aging-brain donors, so these rows are not independent replications. A sign difference between the two cell–gene combinations is not a formal cell interaction test.

Crucially, the previous [regional single-signal coloc pilot](../../../chr5_brain_xqtl/dense_coloc/RESULTS.md) favored **distinct dominant ADHD and expression signals** for both selected cell–gene pairs under its tested priors. Its insomnia comparisons were prior-sensitive; neither finding excludes shared secondary signals or proves a gene mediates either trait. The new exact-variant associations do not change those conclusions or link the reporter assay mechanistically to either gene. rs62369151 is in the ADHD-prioritized chromosome 5 region, distinct from the stronger insomnia–ADHD colocalization region.

The [protocol](PROTOCOL.md), [exact QTL records](exact_qtl_rows.tsv), [portal response](portal_response), [source receipt](receipt.json), [code](lookup.py) and [verification](verification.json) preserve the four-slot analysis. Both source files and tabix indices matched previously recorded publisher MD5s; tabix warned that index timestamps preceded the data file timestamps, but the index hashes matched and lookup succeeded. The source studies use aging-brain tissue; experimental perturbation in relevant human cells, better-matched LD and independent disease data are needed. The failed LAVA gates are unchanged.

Data for this study were prepared, archived and distributed by the National Institute on Aging Alzheimer's Disease Data Storage Site (NIAGADS), University of Pennsylvania (U24-AG041689); accession **NG00184**.

# Reporter-tested variants linked to rs62264764

Among three chromosome 3 variants in the selected Siraj reporter-design table, **rs9815484 and rs1456193 each have r²=0.821 with rs62264764** in the 503-person European reference (bootstrap interval 0.761–0.880 for each). rs73187947 has r²=0.016 (0.006–0.031). The rs1456193 result reproduces the earlier reference check. These are finite reference-panel estimates, not the actual LD in every study population.

The earlier SK-N-SH neural-derived reporter result for rs9815484 had nominal allele-skew p=0.142 and source FDR=0.419. We additionally checked its OL33 rows in four non-neural cell lines:

| Cell line | Reported log₂ allele skew (SE) | Nominal p | Source skew FDR |
| --- | ---: | ---: | ---: |
| GM12878 | −0.118 (0.246) | 0.631 | 0.920 |
| HepG2 | −0.121 (0.218) | 0.579 | 0.856 |
| A549 | +0.128 (0.131) | 0.331 | 0.772 |
| K562 | −0.273 (0.150) | 0.0679 | 0.264 |

None of these five allele-skew rows meets the source's FDR<0.05 criterion. The non-neural rows have both allele-activity Bonferroni values equal to 1; the SK-N-SH row has one lower activity p-value but does not meet the allele-skew threshold. We do not label the marker an author-defined expression-modulating variant based on these rows. The source FDR is not correction for the wider post-hoc Brain6 search. The [rs1456193 follow-up](../rs1456193_other_cells/RESULTS.md) similarly found no source-FDR significant skew in the five checked contexts. rs73187947 appears in the design but has no matched SK-N-SH row in the earlier scan; that absence is not an experimental negative.

These reporter experiments test **other variants**, not rs62264764, in short constructs outside endogenous chromatin. Their nonsignificant results do not establish zero effect, rule out relevant brain or developmental contexts, validate or refute AlphaGenome, assign LSAMP/LINC03051, or establish a shared causal sleep–ADHD locus. Strong LD is useful for follow-up but cannot identify the causal variant; the original local-QC failures are unchanged.

The [pre-result plan](PROTOCOL.md), [reference LD code](analyze.py), [LD table](ld.tsv), [reference receipt](receipt.json), [four-cell query](other_cells/query.py), [source rows](other_cells/source_rows.tsv), [archive receipt](other_cells/receipt.json), [run log](other_cells/run.log), and [independent arithmetic verification](verification.json) preserve the analysis. Source archive members passed ZIP CRC and gzip reading, and exact GRCh37 alleles matched; neither the full archive nor the remotely indexed VCF was whole-file hashed. The LD bootstrap resamples reference individuals; no individual genotypes were saved. Verification checks allele matching and aggregate counts, not study-matched LD, assay calibration or bootstrap regeneration.

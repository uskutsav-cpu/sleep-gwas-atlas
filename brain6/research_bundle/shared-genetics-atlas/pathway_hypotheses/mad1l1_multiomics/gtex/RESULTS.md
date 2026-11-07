# Exact GTEx v8 associations for two linked MAD1L1-region variants

The [GTEx Portal](https://gtexportal.org/api/v2/redoc) v8 single-tissue significant-pair endpoints returned records for the exact GRCh38 A>G alleles of rs28528997 and the linked methylation instrument rs56305291. All four eQTL/sQTL queries succeeded; all returned records, tissue labels and response hashes are in [associations.tsv](associations.tsv), [queries.tsv](queries.tsv) and the saved `responses/` directory. This was a post-selected lookup, not a genome-wide or tissue-enrichment analysis.

| Variant (ALT=G) | Measured eQTL gene | Tissue | Nominal p | ALT normalized effect size |
|---|---|---|---:|---:|
| rs28528997 | AC110781.3 (ENSG00000176349.11) | Nucleus accumbens, basal ganglia | 1.51×10⁻⁶ | −0.492 |
| rs56305291 | AC110781.3 (ENSG00000176349.11) | Nucleus accumbens, basal ganglia | 1.68×10⁻⁵ | −0.454 |
| rs28528997 | MRM2 (ENSG00000122687.17) | Caudate, basal ganglia | 7.37×10⁻⁶ | −0.330 |
| rs56305291 | MRM2 (ENSG00000122687.17) | Caudate, basal ganglia | 3.45×10⁻⁵ | −0.309 |

The queried GTEx label AC110781.3 denotes an annotated transcript antisense to *MAD1L1*; the Ensembl lookup describes its unversioned gene ID as a [novel antisense transcript](https://rest.ensembl.org/lookup/id/ENSG00000176349?content-type=application/json). It is **not the MAD1L1 protein-coding gene**. MRM2 is a separate gene; it is an alternative molecular association in this region, not an established mediator of sleep, ADHD or frailty.

The significant-pair endpoint also reports MAD1L1 splice-phenotype associations in **testis** for both variants and in **sun-exposed skin** for rs28528997; these are splice-cluster observations, not gene-expression results or brain-specific evidence. No MAD1L1 eQTL or brain MAD1L1 sQTL record was returned for these exact alleles by these endpoints. Absence from significant-pair output does not prove absence of a true effect; power, tissue, cell type and variant coverage matter. The rs28528997 query addresses A>G only, while its Ensembl record also lists a T allele.

GTEx significant pairs pass gene-specific discovery thresholds; displayed p-values are nominal and not corrected for the project's post-hoc locus/tissue choices. NES is an ALT-versus-REF regression slope in normalized expression or the stated splice phenotype, not a fold change or disease-risk effect. The two linked variants are measured in overlapping GTEx donors, so their records are not independent molecular replications. A QTL at a selected SNP does not establish that the QTL and either GWAS share a causal variant; dense, multi-signal GWAS–QTL colocalization and experimental validation remain necessary. The blood methylation result from Qian et al. is a separate molecular observation and is not validated by these eQTLs.

[Protocol](PROTOCOL.md) and [query code](query.py) preserve the selection and source requests; the response files contain aggregate association records only. The gene-symbol mapping was checked against [Ensembl MRM2](https://rest.ensembl.org/lookup/id/ENSG00000122687?content-type=application/json) and [MAD1L1](https://rest.ensembl.org/lookup/id/ENSG00000002822?content-type=application/json). See [GTEx methods](https://gtexportal.org/home/methods) for the threshold definitions.

A subsequent [dense regional QTL and coloc analysis](../dense_qtl/coloc/RESULTS.md) compared both expression signals with R13 insomnia and HFRS. Under the stated single-signal priors, H4 was low in all four pairs; coverage and model limitations preclude ruling out secondary shared signals.

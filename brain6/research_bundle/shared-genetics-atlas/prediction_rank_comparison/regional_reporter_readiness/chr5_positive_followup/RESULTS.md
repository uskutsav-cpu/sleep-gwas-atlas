# Reporter-positive chromosome 5 candidate: genetic context

The post-selected reporter row chr5:103904399 G>A maps to **rs12659965**, GRCh38 chr5:104568698 G>A. It has moderate LD with the existing chromosome 5 insomnia–ADHD leads in the 503-person European reference: r²=**0.513** with rs2431108 and **0.511** with rs77960. The respective reference-bootstrap intervals are 0.465–0.563 and 0.462–0.562. LD alone does not establish a common causal variant.

The exact G/A contrast appears in two available phenotype sources:

| Source | Reported marginal association | Interpretation |
| --- | ---: | --- |
| FinnGen R13 registry insomnia | p=2.635×10⁻⁵ | A-allele beta −0.0280, SE 0.00667; nominal source association. |
| Ye UK Biobank Fried frailty score | p=2.6×10⁻¹¹ | Under the supported A1 effect-allele convention, A beta −0.01216 (SE 0.00183); see the [adapter audit](../../../replication/ye_2023/allele_adapter/RESULTS.md). |

The larger p-value is 2.635×10⁻⁵, an unadjusted max-p descriptive conjunction statistic if both marginal tests are valid. It is **not** a corrected discovery: the candidate was selected after the reporter screen and many earlier analyses. The phenotype cohorts' participant and relative overlap have not been independently verified, and FI, FFS and insomnia definitions differ. A same-marker association with both traits does not establish colocalization, independent replication or a gene mechanism.

The candidate was absent by rsID or position from both archived, hash-verified Brain6 insomnia and ADHD regional input files and from the saved FI regional extract. That means **no project ADHD or FI result can be assigned from those files**; it is not proof that the association is absent from the original GWAS or elsewhere. The [Ensembl GRCh38 record](https://rest.ensembl.org/variation/human/rs12659965?content-type=application/json) lists G/A/C at this position, while the queried 1000 Genomes record and the Finnish row have G/A. Exact alleles, reference build and multiallelic representation matter for any future join.

A later [childhood ADHD-symptom lookup](adhd_symptom/RESULTS.md) found a nominal association at the exact G/A contrast (p=0.01575) with a negative A-allele estimate, like the R13 insomnia estimate. This is post-selected related-phenotype context, not independent clinical-ADHD replication or a corrected three-trait discovery.

An [additional-source lookup](other_traits/RESULTS.md) finds nominal MVP insomnia and Finnish hospital-frailty associations at the same contrast. These post-selected results are not a corrected multi-trait discovery, and the Finnish releases overlap.

In the previously checked SK-N-SH reporter row, allele skew passed the source FDR and an allele passed the working activity screen. This is a **candidate proxy**, not a confirmed author-defined emVar or an endogenous chromatin experiment. No AlphaGenome score was generated for this new marker, and no result identifies an effector gene. The original chromosome 5 trait–trait coloc remains model-dependent and its failed LAVA quality gates are unchanged.

The [protocol](PROTOCOL.md), [analysis](analyze.py), [LD table](ld.tsv), [regional source receipt](receipt.json), [full FFS scan](ffs_lookup.py), [FFS receipt](ffs_receipt.json), [Ensembl response](ensembl.json), and [independent arithmetic and hash verification](verification.json) preserve this post-hoc lookup. The FFS scan checked all 8,883,488 rows against the saved source hash; both Brain6 regional inputs matched their frozen hashes. The reference VCF and index were queried remotely without whole-file hashes; bootstrap intervals resample those 503 reference individuals. No genotype records were saved.

A subsequent [full archived source check](full_source_check/RESULTS.md) also found no matching rsID or GRCh37 position in the project raw ADHD and frailty-index files. This strengthens the statement about those copies' coverage but does not establish the marker's biological effect or its coverage in other releases.

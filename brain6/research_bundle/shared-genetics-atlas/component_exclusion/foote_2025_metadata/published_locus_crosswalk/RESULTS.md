# Published frailty-locus comparisons and sleep inputs

[Foote et al. 2025](https://www.nature.com/articles/s41588-025-02269-0) compared previously reported FI and Fried frailty loci with their multivariate factors. Their supplementary tables give the following author-reported classifications for project-relevant leads. “Located within latent-factor locus” and “p<0.001” are the table's own categories; they are not independent replication of a joint sleep–frailty association.

| Earlier frailty lead | Foote supplement | Author-reported category | Factor; reported Q-locus flag |
|---|---|---|---|
| rs1363103, FI | Table 13, row 9 | Located within a latent-factor locus | General factor; Q-locus flag yes |
| rs2396766, FI | Table 13, row 10 | Not located within a latent-factor locus; p<0.001 category yes | Multimorbidity; Q-locus flag no |
| rs8089807, FI, chr18 | Table 13 | Not located within a latent-factor locus; p<0.001 category yes | General factor; Q-locus flag yes |
| rs7703746, FFS | Table 14, row 16 | Located within a latent-factor locus | General factor; Q-locus flag yes |

The table's `Q Locus` classification is a **regional** label and should not be silently equated with the Q p-value of any particular SNP. For example, the separate indexed general-factor record for rs1363103 reports Q p≈0.029. Those describe different quantities; the unfiltered per-SNP and locus construction would be needed to reconcile all details. The authors' “replication” comparison uses their own factor definitions and overlap criteria, and their inputs include UK Biobank GWAS; it is not evidence that independent cohorts replicated both sleep and frailty effects. The table rows are transcribed in `published_locus_rows.tsv` with source row numbers, including the authors' original p and direction fields.

The study's Supplementary Table 1 identifies its tiredness indicator as a Pan-UKB GWAS of field 2080 (408,124 reported) and its insomnia indicator as a UK Biobank binary GWAS (386,988 reported). These are indicators within a 30-trait genomic model, not independent participants counted across traits. This documents concrete phenotype and cohort overlap with the project's UKB-based frailty instruments. The indicators' source details are transcribed in `sleep_input_rows.tsv`; sample counts are those printed for individual source studies, not the number of independent people in the multivariate model.

The reproducible extraction and workbook SHA-256 are in this folder. The public [supplementary workbook](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-025-02269-0/MediaObjects/41588_2025_2269_MOESM4_ESM.xlsx) is read-only source material and is not bundled in full. Its original table columns, not our own significance calculations, determine the classifications above. Regional correspondence of rs8089807 to the chr18 project candidates is proximity only; the table does not prove common causal variants.

## Previously reported chromosome 5 general-factor locus

Supplementary Table 6, row 23, explicitly identifies **rs2431108** as the top lead for general-frailty genomic locus 20, with lead p=1.84136×10⁻²² and locus boundaries chr5:103673480–104092037 (hg19). Its listed independently significant SNPs include rs1363103; that is the paper's LD-based locus terminology, not an independent cohort replication. Thus the general-frailty association at the Brain6 lead was already a published locus, not merely a new association we discovered from the downloaded file. The exact transcribed row is in `published_general_locus_row.tsv`.

The table 13 and 14 rows for rs1363103 and rs7703746 each report `Locus P-value` 7.84619×10⁻⁹, whereas Table 6 reports 1.84136×10⁻²² for the rs2431108 lead. They appear in different comparison tables and are labeled as different table fields; their exact relationship and calculation were not resolved here. Do not substitute either number for the other or call this difference an error without author clarification. The paper's locus boundaries and classification are not proof of a shared causal variant with sleep or ADHD; the factor includes insomnia among its inputs and constituent cohorts may overlap the Brain6 and frailty sources.

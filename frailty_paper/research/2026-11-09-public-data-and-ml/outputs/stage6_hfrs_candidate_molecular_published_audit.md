# Hospital-frailty published molecular results: candidate overlap

We compared gene symbols in the project's targeted molecular inventory (`stage5_variant_gene_tissue_phenotype_integrated.csv`, 98 distinct values) with Supplementary Tables 11–14 of [Mak et al. (2025)](https://www.nature.com/articles/s43587-025-00925-y). This is an exact-name, postselected cross-reference, not a comprehensive search: the inventory and publication tables both apply their own selection, and a missing row is not a null association. The study used a single-causal-variant colocalization model, default priors and 1-Mb windows, and examined FinnGen HFRS with GTEx v8 eQTL/sQTL or UK Biobank Pharma Proteomics Project pQTL. These are HFRS rather than FI or FFS results and do not include an insomnia colocalization test.

| Published comparison | Reported result | Interpretation |
|---|---|---|
| MAD1L1 plasma pQTL versus full HFRS | 14,715 shared SNPs; H3 100.0%, H4 <0.05% | The published model favors distinct signals; precision is limited to the table's rounding. |
| NCAM1 plasma pQTL versus full HFRS | 9,131 shared SNPs; H3 100.0%, H4 <0.05% | Same conditional interpretation. |
| NCAM1 plasma pQTL versus HFRS without dementia | 9,131 shared SNPs; H3 100%, H4 <0.05% | Similar result in a phenotype sensitivity analysis using the same study data. |
| FOXP2 adrenal eQTL versus HFRS without dementia | One shared SNP; H4 0.7547 | Only one shared SNP was analyzed. Distinct two-variant H3 is impossible on this support, so the reported H4 cannot distinguish same-variant sharing from two different causal variants omitted from the intersection. It is not evidence of adrenal mediation. |
| ADARB1 splicing versus HFRS without dementia | Eight tissue entries with only 1–2 shared SNPs; reported H4 often >0.98 | High H4 on this limited support is fragile; a one-SNP comparison cannot model distinct causal variants, and even two provide little coverage. |

For MAD1L1 and NCAM1 the posterior ranking is conditional on the paper's model, priors, allele harmonization and variant coverage, and does not prove that no molecular mechanism exists. Both protein and HFRS sources share populations with other project analyses, so these rows are not independent validation of a candidate or score. No gene was matched in the full-HFRS GTEx eQTL/sQTL table; this is a property of the selected rows and exact name matching only.

Reproducibility: `work/stage6_hfrs_candidate_pqtl_audit.py` reads the previously saved workbook `work/hfrs_2025_tables.xlsx` and produces `outputs/stage6_hfrs_candidate_pqtl_audit.json`. The input workbook's recorded SHA-256 is `123c340cf063e2ef618b8ef543dedea6b7c8ef69008488f01027eedbffc04625` (see `rs12536335_finngen_hfrs_followup.md`). The 10 dementia-excluded GTEx rows also contain one SSR4P1 row with four shared SNPs, subject to the same limited-support concern.

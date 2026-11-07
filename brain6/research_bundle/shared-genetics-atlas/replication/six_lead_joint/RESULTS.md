# Six original insomnia–ADHD leads in two follow-up sources

Of the six original Brain6 insomnia–ADHD regional leads, **five matched both the FinnGen R13 insomnia file and the childhood ADHD-symptom GWAS**. The chromosome 5 lead rs2431108 showed the clearest cross-source consistency: its insomnia p=0.000124 and symptom p=0.00600, with both effects aligned in the same directions as their respective Brain6 discovery effects. This is a promising **related-phenotype follow-up**, not yet confirmatory independent replication.

| Original lead | FinnGen insomnia p | Childhood symptoms p | max(p_sleep, p_symptom) | max-p × 6 |
|---|---:|---:|---:|---:|
| rs2582894, chr11 | 0.0158 | 0.195 | 0.195 | 1.000 |
| rs7160630, chr14 | 0.0540 | 0.806 | 0.806 | 1.000 |
| rs62264767, chr3 | 0.00592 | 0.478 | 0.478 | 1.000 |
| **rs2431108, chr5** | **0.000124** | **0.00600** | **0.00600** | **0.0360** |
| rs6452785, chr5 | No matching FinnGen record | 0.0632 | Not estimable | Not estimable |
| rs1476535, chr7 | 0.00417 | 0.660 | 0.660 | 1.000 |

The six rows are the actual original region leads, not six independently selected follow-up hits. For rs6452785, the separate scan returned no FinnGen row at the lifted coordinate; this is missing evidence, not a negative association. It remains in the denominator of six. rs62264767 is near the AlphaGenome-prioritized rs62264764 but is a different variant. The ADHD-symptom association at rs62264767 is weak, so this analysis does not substantiate an ADHD-symptom association at that chr3 lead.

For rs2431108, aligned to the C allele, the R13 insomnia beta was +0.02811 (SE 0.00732), and the childhood symptom beta was +0.00809 (SE 0.00294). Sleep is modeled as log odds; symptoms use a quantitative score, so their effect magnitudes are not directly comparable. In 503 European reference participants, rs2431108 and the previously reviewed rs77960 have r²≈0.995; the two observations should not be counted as independent loci.

## Interpretation of the six-lead benchmark

If both marginal p-values are valid, max-p is an intersection–union test of association with both phenotypes and does not require independence between the two follow-up statistics. Multiplying by six is the Bonferroni benchmark for this specified six-lead set. **It is not a confirmatory adjusted replication result here:** the chromosome 5 FinnGen result and related rs77960 follow-up were previously inspected, the original discovery family was larger, and the retrospective selection and earlier analyses are not captured by a factor of six. Direction checks are descriptive; weak associations do not firmly establish an effect direction.

The [childhood ADHD-symptom study](https://www.nature.com/articles/s41588-025-02295-y) reports that its symptom sample and the Demontis diagnostic ADHD discovery GWAS did not overlap. That is author-reported participant nonoverlap, not an independent audit of identity, relatedness or overlap with all sleep cohorts. The symptom measure is not a clinical diagnosis. The [FinnGen insomnia endpoint](https://r13.risteys.finngen.fi/endpoints/F5_INSOMNIA) is registry-defined and differs from the original questionnaire phenotype; exact discovery–follow-up participant and relative overlap remains unverified. The two follow-up sources do not establish one causal variant or an effector gene. A prospectively frozen full-family analysis with source-verified cohorts and suitable LD remains needed; the failed LAVA gate is unchanged.

## Reproducibility

[Protocol](PROTOCOL.md), [six-lead table](six_leads.tsv), [full source receipt](receipt.json), [missing-coordinate check](missing_coordinate_check.json), [code](analyze.py) and [independent arithmetic check](verification.json) record the selection, alignments, matching and data lineage. Original Brain6 regional inputs matched frozen SHA-256 hashes. The entire 21,327,010-row FinnGen gzip was read and matched the previously recorded source SHA-256; the symptom GWAS was queried through the [indexed GWAS Catalog file](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90568001-GCST90569000/GCST90568440/harmonised/) without full-file or index hashes. The verification checks family membership, symptom orientation and max-p arithmetic; it does not validate sample independence or selection adjustment. No AlphaGenome predictions were used in the tests.

A subsequent [MVP insomnia check of all six leads](../mvp_six_leads/RESULTS.md) found rs2431108 p=0.0895 despite the same effect direction, weakening cross-sleep-source consistency. The different phenotype and power mean this is not proof of no effect.

An earlier [falls GWAS follow-up](../chr5_falls_2020/RESULTS.md) reports the same chromosome 5 variant in UK Biobank, Rotterdam and B-PROOF. Its mixed cohort-specific results concern falls, not independent replication of the insomnia–ADHD pair or the FI/FFS/HFRS phenotypes.

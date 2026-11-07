# MVP insomnia follow-up of six Brain6 insomnia–ADHD leads

A third sleep-source check **weakens the cross-source support for chromosome 5 rs2431108**. Its Million Veteran Program (MVP) insomnia p-value is **0.0895**, with the same effect direction as Brain6's original insomnia GWAS. This is weaker evidence than its FinnGen registry-insomnia p=0.000124. A weak association does not establish that there is no effect; phenotype, population, model, precision and sample differences can all matter.

| Original lead | MVP insomnia p | Childhood ADHD-symptom p | max-p | max-p × 6 |
|---|---:|---:|---:|---:|
| rs2582894, chr11 | 0.000933 | 0.195 | 0.195 | 1.000 |
| rs7160630, chr14 | 0.000791 | 0.806 | 0.806 | 1.000 |
| rs62264767, chr3 | 0.0000499 | 0.478 | 0.478 | 1.000 |
| rs2431108, chr5 | 0.0895 | 0.00600 | 0.0895 | 0.537 |
| rs6452785, chr5 | 0.0000310 | 0.0632 | 0.0632 | 0.379 |
| rs1476535, chr7 | 0.000000120 | 0.660 | 0.660 | 1.000 |

All six original leads matched exactly by GRCh38 position, rsID and allele pair in both sources. None crosses the displayed 0.05 numerical benchmark. This is a **post-selected descriptive comparison**: related candidates and other cohorts were already inspected, and the original discovery family was larger than six. The max-p construction is an intersection–union test when both marginal p-values are valid and does not require independence between the two follow-up traits, but multiplying by six here is not a confirmatory correction for prior selection. A non-crossing does not exclude a shared biological effect. Weak childhood-symptom evidence, rather than the MVP insomnia p-value, limits the max-p for five leads; rs2431108 is limited by its MVP insomnia p-value.

MVP's source reports a rounded T-allele odds ratio of 0.9886 for rs2431108, corresponding to an aligned C-allele log-odds estimate of approximately +0.0115. The supplied SE field is missing; we use the source p-value directly, not a p-value reconstructed from the rounded OR. The childhood symptom effect on the C allele is +0.00809, with p=0.00600. Those signs match their respective original Brain6 discovery effects, but the outcomes and effect scales differ. The prior FinnGen result and this MVP result are separate cohort-level evidence, **not two independent disease loci** or proof that the same causal variant drives either phenotype.

The source [GCST90475826](https://www.ebi.ac.uk/gwas/studies/GCST90475826) is MVP European-ancestry insomnia; the study uses clinical-record ascertainment, while the original sleep GWAS and the FinnGen registry endpoint differ in measurement and populations. The childhood study reports no participant overlap with the 2023 diagnostic ADHD GWAS, but it measures symptoms rather than a clinical diagnosis. Exact participant or relative overlap among the original, MVP, FinnGen and symptom cohorts has not been independently audited. The chr5 finding is therefore best treated as **inconsistent strength across sleep follow-ups**, not validated two-trait replication or a negative result. No local-correlation gate, causal gene or AlphaGenome prediction is validated by this analysis.

The [protocol](PROTOCOL.md), [six-lead values](six_leads.tsv), [source and query receipt](receipt.json), [analysis code](analyze.py), and [independent arithmetic check](verification.json) preserve this post-hoc screen. The complete 19,703,815-row MVP gzip matched the previously recorded compressed-file SHA-256. Discovery regional input hashes and exact allele/coordinate matching were checked. Source p-values and rounded ORs were not used to derive each other. Arithmetic and file-identity checks do not validate cohort independence or correct for previous selection. The [earlier FinnGen comparison](../six_lead_joint/RESULTS.md) remains useful context with its own limitations.

A subsequent [clinical-insomnia effect-difference sensitivity](../mvp_finngen_clinical/sleep_effect_sensitivity/RESULTS.md) compares the two estimates directly. Unknown source SE and cohort dependence prevent a calibrated heterogeneity conclusion.

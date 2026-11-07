# Cross-source sleep and childhood ADHD-symptom follow-up

At **rs77960**, the existing FinnGen R13 insomnia data and the separate childhood ADHD-symptom GWAS both show associations in the same allele directions as their respective Brain6 discovery effects. Their nominal p-values are **0.000347** for registry insomnia and **0.00606** for childhood symptoms. This is a promising cross-source consistency result for a previously published chromosome 5 region, but it is **not yet confirmed independent two-trait locus replication**.

| Variant, aligned allele | FinnGen insomnia beta (SE); p | Childhood symptom beta; p | max(p_sleep, p_symptom) | max-p × 3 |
|---|---|---|---:|---:|
| rs62264764, A | −0.02248 (0.00813); 0.00569 | −0.00291; 0.437 | 0.437 | 1.000 |
| rs77960, A | +0.02627 (0.00734); 0.000347 | +0.00808; 0.00606 | 0.00606 | 0.0182 |
| rs1476535, C | −0.01914 (0.00668); 0.00417 | −0.00119; 0.660 | 0.660 | 1.000 |

Both trait signs agree with their corresponding Brain6 discovery effects for all three markers, although the two weak symptom associations are not persuasive evidence of effect direction. The sleep betas are log-odds estimates; the symptom betas are on a quantitative scale and should not be compared numerically. The `max-p` statistic is an intersection–union test if each marginal p-value is valid, without requiring independence between the two follow-up traits. Multiplication by three counts the three variants in this retrospective comparison. Because rs77960 and its FinnGen result were already seen before this integration, and the broader candidate family and earlier searches were larger, **0.0182 is a descriptive benchmark, not a confirmatory family-wise adjusted p-value**.

## What independence would require

The [childhood symptom study](https://www.nature.com/articles/s41588-025-02295-y) reports no sample overlap with the diagnostic ADHD GWAS used for Brain6, giving author-reported nonoverlap on that side. The ADHD follow-up is a childhood symptom measure rather than an independent diagnosis GWAS. [FinnGen R13 F5_INSOMNIA](https://r13.risteys.finngen.fi/endpoints/F5_INSOMNIA) is a registry-defined insomnia endpoint, not the original Brain6 questionnaire phenotype. FinnGen and the original sleep cohort come from distinct described sources, but their exact participant intersection and relatedness have not been independently checked; nor has overlap between the follow-up traits. The FinnGen values were already used in earlier frailty analyses in this atlas, so this integration is not a fresh blinded validation. It does not establish that the same causal variant drives both traits, identify a gene, repair failed LAVA quality control, or prove a mechanism. A prespecified full-family evaluation with source-verified independent sleep and diagnostic ADHD cohorts and appropriate LD remains needed.

## Source and verification

The [protocol](PROTOCOL.md) explicitly records that this comparison was post hoc. The [three-row comparison](comparison.tsv), [code](compare.py), [input hashes](receipt.json) and [verification](verification.json) preserve the calculations. Exact rsID, GRCh38 coordinate and allele pair match in the previously saved R13 regional extractions, whose hashes match their extraction receipts. The full raw R13 file was previously matched to the official object by MD5; the childhood GWAS was accessed by bounded indexed lookup, without a full-file hash. The source-value and max-p arithmetic checks do not verify cohort independence, account for selection or calibrate the GWAS. No AlphaGenome score was used.

A subsequent [six-original-lead comparison](../six_lead_joint/RESULTS.md) includes the closely linked chromosome 5 original lead rs2431108 and five other original leads. The six-lead benchmark remains retrospective and is not a confirmatory replication correction.

A later [MVP insomnia follow-up](../mvp_six_leads/RESULTS.md) found a weak rs2431108 association (p=0.0895) for the closely linked original chr5 lead. This constrains claims of sleep-side consistency and does not establish absence of an effect.

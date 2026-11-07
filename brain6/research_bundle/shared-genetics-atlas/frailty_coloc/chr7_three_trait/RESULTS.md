# FOXP2 region: FI, insomnia and ADHD common-support pilot

On 1,055 matched chromosome 7 variants, the default exploratory HyPrColoc model returned an **FI–insomnia cluster**, not a three-trait cluster. The cluster posterior was **0.8822** under the stated model and independence assumption. More stringent settings weakened or removed the cluster. This is not evidence that ADHD has no association or no secondary shared signal.

| Setting | Returned cluster | Model-conditional posterior | Candidate variant and conditional contribution |
|---|---|---:|---|
| Conditional prior 0.02; thresholds 0.5 | FI, insomnia | 0.8822 | rs12671330; 0.0811 |
| Conditional prior 0.01; thresholds 0.5 | FI, insomnia | 0.7870 | rs12671330; 0.0811 |
| Conditional prior 0.005; thresholds 0.5 | FI, insomnia | 0.6428 | rs12671330; 0.0811 |
| Conditional prior 0.02; thresholds 0.9 | None | Not reported | None |

The conditional contribution is the fraction of that cluster posterior attributed to the candidate under the model, **not** an unconditional probability of variant causality. The association prior is 10⁻⁴; the other prior is the assumed chance of an additional trait sharing a variant. Cluster posteriors are not p-values or multiple-testing-corrected discoveries.

## Pairwise context on the same 1,055 variants

| Pair | H3, distinct, p12=10⁻⁵ | H4, shared, p12=10⁻⁵ | H4, p12=10⁻⁶ |
|---|---:|---:|---:|
| FI–original Brain6 insomnia | 0.0442 | 0.9528 | 0.6688 |
| FI–ADHD | 0.9557 | 0.0441 | 0.00460 |
| Original Brain6 insomnia–ADHD | 0.1450 | 0.8454 | 0.3535 |

These single-signal ABF posteriors are conditional on priors, selected support and each pair's model. They cannot be multiplied into a three-trait posterior. In particular, the FI–ADHD comparison favors distinct *dominant* signals in this model, not proof of distinct causal variants or absence of a secondary shared signal. The earlier FI–FinnGen registry-insomnia model used a different sleep source and 1,212 markers and had H4=0.5873 at p12=10⁻⁵; cohort, phenotype, model and support all differ, so the comparison does not identify the cause of the difference.

## Reference linkage and phenotype overlap

The cluster marker rs12671330 has r²=0.110 with Brain6 marker rs1476535 and r²=0.106 with frailty candidate rs12536335 in 503 European 1000 Genomes reference participants. It has r²=0.714 and 0.747 with the two prior FI–FinnGen insomnia conditional-H4 leaders rs10953765 and rs1456031. These are post-hoc reference-panel measurements, not study-matched LD, a conditional test, or evidence of independent causal signals. The conditional leader is not automatically an effector gene; FOXP2 is a regional label.

The original Brain6 insomnia GWAS and FI share UK Biobank insomnia measurement content, and FI also contains tiredness. FI and other cohorts may overlap; independence of the three association statistics is **unverified**, while HyPrColoc was run assuming independence. Thus the returned FI–insomnia cluster could partly reflect measurement and sample overlap and does not show a sleep-independent frailty mechanism. FI meta-analysis cohort conventions, including the unresolved TwinGene allele-adapter question, ancestry differences, rounded FI statistics, source GWAS frequency omissions, and the previously selected FinnGen-filtered support remain limitations. The analysis has no FFS or HFRS and is not a sleep-item-excluded FI, independent replication, a local genetic correlation, a validated multi-signal fit, or a rescue of failed LAVA QC. The models assume at most one causal variant per trait; they can misrepresent regions containing multiple signals. AlphaGenome predictions were not model inputs.

[Protocol](PROTOCOL.md), [matching receipt](matching_receipt.json), [pairwise posteriors](pairwise.tsv), [cluster outputs](clusters.tsv), [reference LD table](pairwise_ld.tsv), [model code](run.R) and [verification](verification.json) document the analysis. The package's documented example reproduced its expected candidate markers; cluster SNP scores normalized and matched the input order. An independent implementation reproduced all six pairwise H0–H4 vectors within 1.5×10⁻¹⁵. These are arithmetic and software checks, not calibration or biological validation. The restricted source-derived input stays on the SSD; the PGC ADHD source's noncommercial and redistribution terms apply. The reference-panel bootstrap intervals in `pairwise_ld.tsv` describe that panel only.

A [precision stress test](precision_sensitivity/RESULTS.md) kept the default FI–insomnia cluster in all 200 illustrative FI rounding perturbations, but its top variant changed in 81 of them. These perturbations do not establish causal-marker probabilities or scientific validation.

A subsequent [illustrative study-correlation stress test](overlap_sensitivity/RESULTS.md) changed the returned cluster from FI–insomnia to insomnia–ADHD for assumed positive FI–insomnia correlations of 0.1 or higher. The real correlation is unknown, so the result demonstrates model sensitivity rather than identifying the true cluster.

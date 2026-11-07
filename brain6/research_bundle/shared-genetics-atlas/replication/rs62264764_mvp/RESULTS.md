# rs62264764: MVP insomnia and the nearby original lead

The exact AlphaGenome-prioritized variant **rs62264764** is associated with MVP clinical-record insomnia (source p=**4.918×10⁻⁵**), with its A-allele effect in the same direction as in the Brain6 discovery insomnia GWAS. The earlier FinnGen registry-insomnia lookup also had the same direction (p=0.00569). These post-selected associations across sleep sources are useful consistency evidence, but exact cohort independence and prior-selection correction remain unresolved.

| Marker | MVP insomnia p | FinnGen insomnia p | Childhood ADHD-symptom p |
|---|---:|---:|---:|
| rs62264764 | 0.00004918 | 0.00569 | 0.437 |
| Nearby original lead rs62264767 | 0.00004993 | 0.00592 | 0.478 |

These are **marginal** associations in different datasets and phenotypes, not a shared-causal-variant test or independent two-trait replication. Weak childhood-symptom evidence does not establish the absence of an ADHD effect. The source MVP odds ratio is rounded and its SE field is unavailable, so the source p-value was used directly; aligned A-allele log(OR) is approximately −0.0315. The symptoms GWAS authors report no participant overlap with the diagnostic ADHD source, but exact overlap and relatedness among all sleep and follow-up cohorts have not been independently verified. The marker and sources were selected after earlier analyses, and these p-values are not corrected for the wider research search.

A later [larger ADSP European-labeled reference check](../../variant_context/rs62264764_brain_xqtl/ldz_readiness/adsp_public_reference/RESULTS.md) reports r² = 0.99919 for the same pair. It strengthens the reference-panel tagging observation but is not disease replication or a study-matched full LD matrix.

## Why the nearby lead matters

In the 1000 Genomes Phase 3 European reference, **rs62264764 A and rs62264767 C were identical across all 1,006 phased haplotypes from 503 people** (r²=1.0 in this panel). This is a finite reference observation: phase, panel sampling and the study populations' LD can differ. It does not prove that either variant is causal or that they are indistinguishable in every cohort. Their similar GWAS effects are expected for strongly tagged markers but do not identify an effector gene.

The saved AlphaGenome adult frontal-cortex DNase-track predictions differ: the 501-bp ALT-versus-REF log₂ ratio is **−0.2734** at rs62264764 and **+0.0130** at rs62264767. These are predictions in windows centered at different variants, not measured effects or disease-risk probabilities. The contrast suggests a testable allele-specific accessibility experiment at both variants, but it cannot choose the causal variant, assign an LSAMP effect, or establish brain-cell mediation. It is especially important to test the two variants separately in a biologically relevant cell type and compare predicted and measured accessibility and gene expression. AlphaGenome outputs remain subject to the user's authorized non-commercial research terms and were not used to train a model here.

A [five-group reference LD follow-up](cross_population_ld/RESULTS.md) found r² from 0.987 to 1.000 in the sampled groups. This is not a disease association, independent replication or evidence of identical LD in the study cohorts.

The [protocol](PROTOCOL.md), [MVP result](mvp_result.tsv), [exact source receipt](query_receipt.json), [reference LD table](pairwise_ld.tsv), [prediction comparison](alpha_comparison.tsv), [verification](verification.json) and code record the analysis. The complete MVP gzip matched its previously saved SHA-256; the remotely queried 1000 Genomes VCF and index were not downloaded or whole-file hashed. The 2,000-replicate bootstrap resamples individuals in the reference panel and is not an independent cohort or cross-ancestry validation. These checks do not resolve the chromosome 3 LD–Z diagnostics, establish secondary colocalization, or repair failed LAVA quality control.

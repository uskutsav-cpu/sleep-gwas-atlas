# Connecting Brain6 and frailty candidates: reference-panel LD

Two sleep–frailty candidate markers fall in regions also selected by Brain6 for insomnia–ADHD. To determine whether proximity represents substantial linkage in a reference panel, we queried the indexed [1000 Genomes Phase 3 GRCh37 VCF hosted by UCSC](https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/) for each exact rsID, position and allele pair. We used 503 European-labeled reference participants and 1,006 phased haplotypes for each pair; all listed two-allele contrasts were present in the queried VCF records. Correlations below are between the **named allele indicators**, not between trait effects. Confidence intervals are percentile bootstrap intervals from 2,000 resamples of people (both haplotypes together), intended to describe this panel's sampling variation only.

| Brain6 marker (counted allele) | Frailty candidate (counted allele) | Reference r | r² | Bootstrap r² 95% interval |
|---|---|---:|---:|---:|
| rs62264764 A | rs817509 G | −0.839 | 0.703 | 0.628–0.777 |
| rs1476535 C | rs12536335 G | +0.978 | 0.957 | 0.931–0.980 |
| rs9969232 G | rs12536335 G | +0.785 | 0.617 | 0.564–0.674 |

The chr3 and chr7 positional overlap was identified before this LD check; selecting the pairs post hoc and related reference individuals limit inference. A high r² indicates that the markers tag similar variation in this European reference; it does **not** establish that they are the same causal variant, that two GWAS share a causal signal, or that a gene or biological pathway mediates the phenotypes. The panel does not contain the source GWAS outcomes, may not represent the actual GWAS populations, and phased haplotypes can contain error. The bootstrap assumes the sampled people are independent and is not a cross-ancestry or cohort-replication interval.

## Chromosome 3: LSAMP region

Brain6 uses GRCh37 rs62264764 G>A at 117,639,575 (GRCh38 chr3:117920728 G>A). Its A allele has negative estimated insomnia and ADHD effects in the source GWAS. rs817509 is at GRCh37 position 117,750,861 under the Song candidate build assumption; its G allele has positive estimated FI, FFS, UKB insomnia and HFRS effects in the project's post-selected source lookup. The alleles have negative reference r, so these reported signs are **directionally compatible as a hypothesis**. The source scales differ, some cohorts overlap, and the rs817509 FFS convention is inferred. This does not demonstrate one shared three-phenotype causal effect. rs62264764's GTEx adrenal LSAMP eQTL and AlphaGenome frontal-cortex prediction are additional observations of different types; neither validates the genetic bridge or identifies an effector tissue.

## Chromosome 7: FOXP2 region

Brain6's rs1476535 C (GRCh37 114,071,035) is in strong positive reference LD with frailty candidate rs12536335 G (114,043,159). Both counted alleles have negative estimated insomnia effects in the respective source tables; rs1476535 C is negative for ADHD and rs12536335 G is negative for the FI meta-analysis, FFS and HFRS in the project's descriptive review. This is consistent with a regional follow-up, not an established shared causal variant. The rs12536335 FI cohort estimates disagree in direction under a TwinGene allele adapter that is supported by a broad reconstruction audit but conflicts with the source README and remains unconfirmed by the authors. In addition, an apparent match across sources that reuse participants is not independent replication. rs9969232 provides another nearby ADHD-credible-set marker with r²≈0.617 to rs12536335; it should not be counted as a separate independent locus.

The user-provided project candidate rsIDs were selected after earlier analyses, and Song's source build is not independently authenticated; the exact liftover and study alleles still require source confirmation. rs817509 has additional alleles in broader Ensembl reference records even though the queried 1000 Genomes record is biallelic. GWAS-specific LD, fine-mapping with valid LD, multi-signal colocalization across *all* traits, and independent cohorts are needed to substantiate either connection. A simple pairwise LD comparison does not substitute for conditional association or colocalization.

## Reproducibility

`calculate_ld.py` queries the cited VCF with pysam, uses the 1000 Genomes Phase 3 sample panel at the path in the script, verifies rsID/position/alleles and phased genotype calls, and calculates r from allele frequencies and joint frequency. `reference_receipt.json` records allele representations, query URLs, panel hash and counts; `pairwise_ld.tsv` contains estimates and bootstrap intervals. No participant-level genotypes were saved. The remote VCF is not downloaded in full or byte-hashed; future indexed responses can change. The companion `six_candidates_six_trait_direction_matrix.csv` is a copied, post-selected source table, not an independently recomputed GWAS.

## Exact frailty markers in the Brain6 pilot and GTEx

Both frailty candidates were present in the original 425-variant pilot; no additional AlphaGenome requests were needed. Their model outputs are saved in `cross_project_brain6_variants.tsv`:

| Variant | Insomnia / ADHD credible-set membership | Insomnia / ADHD conditional PIP | Adult local DNase ALT-versus-REF log₂ ratio | AVI PHRED |
|---|---|---:|---:|---:|
| rs817509 | yes / yes | 0.00474 / 0.000951 | +0.0267 | 0.119 |
| rs62264764 | yes / yes | 0.0402 / 0.0312 | −0.2734 | 20.486 |
| rs12536335 | yes / no | 0.0338 / 0.00286 | +0.0194 | 0.291 |
| rs1476535 | yes / yes | 0.0732 / 0.0163 | −0.0207 | 8.274 |

PIPs are conditional, single-signal values and do not use the frailty data. The four variables are selected on prior information, and differences in AlphaGenome scores are neither statistical tests nor evidence that the higher-scoring variant is causal. Model outputs and AVI reflect different properties from PIP, and ALT/REF prediction direction must not be conflated with the source A1 effect direction.

A direct GTEx v8 significant-pair lookup for exact GRCh38 rs817509 A>G returns adrenal eQTL records for LSAMP (nominal p=2.17×10⁻⁸, ALT NES −0.646) and RP11-384F7.2 (p=2.25×10⁻⁶, NES −0.543). These correspond to the **same GTEx donors and linked region** as the rs62264764 adrenal associations, so they are not two independent molecular replications. rs817509 G has negative reference LD with rs62264764 A; their adrenal LSAMP NES signs are correspondingly opposite. This is consistent with tagging the same regional association, but cannot establish the causal variant, molecular mediator, or brain relevance.

Exact rs12536335 A>G has GTEx FOXP2 eQTL records in pituitary, adrenal and subcutaneous adipose; its adrenal p=5.51×10⁻⁵ and ALT NES −0.366. The significant-pair queries returned no brain records for these two candidate variants. Four dynamic exact-variant tests of rs817509–LSAMP and rs12536335–FOXP2 in cortex and frontal cortex gave nominal p-values between 0.20 and 0.74. They do not exclude effects in other cells, tissues, ages or underpowered samples. The dynamic output summarizes counts and fitted values but deliberately does not save the individual-level arrays returned by the public service. `gtex_query_receipt.tsv`, `gtex_associations.tsv`, `gtex_dynamic_brain_summary.tsv` and saved static responses give exact IDs, tissue context and source hashes. GTEx variant–gene p-values are nominal within significant-pair selection; these follow-up selections are post hoc, and neither eQTL nor signed LD is a GWAS–QTL colocalization test.

## Sensitivity across the highest-PIP candidates

To avoid relying on a single chosen pair, we also queried each region's four highest insomnia PIPs, four highest ADHD PIPs, largest absolute DNase-change variant, highest AVI variant and frailty anchor, deduplicating overlaps. All **19 variant–anchor rows**, including the two anchor self-comparisons, matched the specified rsID, GRCh37 coordinate and allele pair in the reference VCF, and each used 503 people. The full `credible_representative_ld.tsv` records signed r, r², selected reasons, PIPs and credible-set membership; no additional significance tests or causal probabilities were computed.

In chromosome 3, rs817509 has r²≈0.58 with rs1456193 and rs9815484, the two highest insomnia-PIP representatives; r²≈0.70 with rs62264767 and rs62264764; and r²≈0.40 with the highest ADHD-PIP representative rs73163723. In chromosome 7, rs12536335 has r²≈0.96 with top insomnia-PIP rs1476535 but only ≈0.11–0.13 with the next three top insomnia-PIP representatives. Its reference r² with the three highest ADHD-PIP representatives is ≈0.62–0.63. The representative choices are trait-specific PIP rankings within the **selected credible-set union**, not full-locus causal probabilities. For three already examined variant pairs, this separate script reproduced the earlier signed r and r² to 10⁻¹².

The differing LD values show why one linked pair does not resolve the whole region: multiple plausible variants remain, and the fine-mapping here assumes one causal variant per trait. Pairwise reference LD cannot distinguish separate causal signals or replace conditional analysis; the European panel, source allele and cohort caveats above continue to apply. `credible_representative_ld.py` and `credible_representative_receipt.json` document selection, the exact reference alleles, input hashes and query status. No individual genotype data were saved.

## Chromosome 7: conditional coloc leaders versus Brain6 anchors

The two highest conditional-H4 variants in the separate FI–FinnGen insomnia single-signal coloc are rs10953765 and rs1456031 (conditional SNP probabilities 0.2154 and 0.2151, **given H4**). Their high reference LD (r²=0.964; bootstrap interval 0.941–0.984) indicates that they tag very similar variation in this panel. Their counted alleles A and T have negative signed r (−0.982), consistent with opposite signs for the two source effects; neither variant is independent evidence of a second shared signal.

| Pair | Reference r² | Bootstrap 95% interval |
|---|---:|---:|
| rs10953765–rs1476535 (Brain6) | 0.140 | 0.101–0.186 |
| rs10953765–rs12536335 (frailty candidate) | 0.134 | 0.096–0.178 |
| rs1456031–rs1476535 (Brain6) | 0.139 | 0.099–0.183 |
| rs1456031–rs12536335 (frailty candidate) | 0.132 | 0.094–0.176 |

These post-selected reference comparisons suggest that the conditional-H4 leaders tag a different pattern of regional variation from the previously highlighted anchors. They **do not prove independent causal signals or disprove sharing**: pairwise LD cannot perform conditional association, the coloc assumes one signal per trait and its overall H4 is only 0.587 under the default prior and 0.125 at the stricter p12. Source samples differ from the reference and any sample overlap is unresolved. The two leading variants' counted-allele frequencies differ between the reference and source summaries by about 0.075–0.087 for rs10953765 and 0.074–0.081 for rs1456031; the two source summaries are relatively concordant with each other. The modeled files matched exact position and alleles and passed the pre-existing 0.15 source-frequency-difference filter, but these reference differences reinforce the need for study-matched LD and a source-level allele/frequency audit rather than assuming this reference represents each cohort.

The [code](chr7_coloc_top_ld.py), [pairwise table](chr7_coloc_top_ld.tsv), and [receipt](chr7_coloc_top_receipt.json) document the posterior selection, reference inputs, allele coding and sampling. Bootstrap intervals describe reference-panel sampling under the individual-resampling assumption, not source-cohort uncertainty or replication. The original coloc remains exploratory; no local-correlation QC gate has been repaired.

## Reference population sensitivity for the chr7 signals

The apparent source–reference frequency differences above are substantially smaller in relevant 1000 Genomes subgroups. For rs10953765 A, FinnGen reports 0.5618 and the FIN reference is 0.5707; the FI summary reports 0.5501 and the GBR reference is 0.5440. For rs1456031 T, the corresponding frequencies are 0.4404 and 0.4242 for FinnGen/FIN, and 0.4481 and 0.4560 for FI/GBR. The source FI is a meta-analysis, so GBR is a sensitivity comparator rather than a direct representation of the entire FI study; FIN is likewise not the actual FinnGen sample. R13 AF is aligned to the FI effect allele in the saved coloc matching table.

| Reference group | People | rs10953765 A frequency | rs1456031 T frequency | r² between leaders | r² rs10953765–rs1476535 |
|---|---:|---:|---:|---:|---:|
| EUR pooled | 503 | 0.475 | 0.522 | 0.964 | 0.140 |
| FIN | 99 | 0.571 | 0.424 | 0.980 | 0.079 |
| GBR | 91 | 0.544 | 0.456 | 1.000 | 0.173 |
| CEU | 99 | 0.510 | 0.490 | 0.960 | 0.199 |
| IBS | 107 | 0.383 | 0.612 | 0.904 | 0.147 |
| TSI | 107 | 0.388 | 0.607 | 0.981 | 0.092 |

The pooled EUR value combines groups with markedly different allele frequencies here; that is a plausible contributor to its difference from the sources, not proof of the cause or proof of correct source allele coding. The comparison is post hoc and the 91–107-person subgroup estimates are imprecise; a sample r² of 1.000 does not establish perfect population LD. Across these five reference subgroups the leaders are strongly linked to each other and less strongly linked to the Brain6 anchor, but this descriptive pattern is not a conditional association test, population-matched fine-mapping, or independent replication. It does not resolve source allele conventions, overlapping participants, or the failed local-correlation QC gates.

The [frequency table](chr7_ancestry_frequencies.tsv), [all five pairwise subgroup LD results](chr7_ancestry_ld.tsv), [script](chr7_ancestry_sensitivity.py) and [receipt](chr7_ancestry_receipt.json) document the reference sample counts, counted alleles and panel/source hashes; no per-person genotypes are saved. The panel is a public 1000 Genomes Phase 3 reference, and the subgroup samples are not equivalent to the GWAS cohorts.

# Chromosome 5: GTEx expression follow-up for RP11-6N13.1

The previously reported insomnia–ADHD region contains a **testis expression signal** for RP11-6N13.1 (GTEx v8 gene ID ENSG00000251574.6) that shows high shared-signal support with each GWAS in an exploratory single-signal model. This is additional molecular context for the chromosome 5 region, **not evidence that the gene acts in brain or causes either phenotype**. Both analyses use the same GTEx testis signal, and the GWAS inputs are the existing Brain6 discovery data; these are not independent replications or a three-trait colocalization.

| Expression tissue and pair | Matched variants | H3, distinct | H4, shared | H4 with stricter shared prior |
|---|---:|---:|---:|---:|
| Testis expression–insomnia | 1,673 | 0.0433 | 0.9566 | 0.6881 |
| Testis expression–ADHD | 1,633 | 0.0175 | 0.9825 | 0.8490 |
| Frontal cortex expression–insomnia | 0 | Not estimable | Not estimable | Not estimable |
| Frontal cortex expression–ADHD | 0 | Not estimable | Not estimable | Not estimable |

The default uses p1=p2=10⁻⁴ and p12=10⁻⁵; the final column uses p12=10⁻⁶. Posterior probabilities are conditional on the specified model, variants and priors; they are not probabilities that the gene mediates sleep or ADHD. Changing both effect-prior SDs to 0.1 or 0.3 gives H4 ranges 0.9567–0.9567 for insomnia and 0.9826–0.9826 for ADHD at p12=10⁻⁵, and 0.6882–0.6884 and 0.8497–0.8499, respectively, at p12=10⁻⁶. The full table retains all digits. The frontal-cortex source returned no nominal rows for this gene in the selected interval; that is **missing or unavailable evidence**, not a negative brain association or proof the gene is unexpressed.

## How this relates to existing results

The testis single-variant associations at rs2431108 and rs77960 agree with the previously saved GTEx significant-pair p-values and normalized slopes. Their respective nominal expression p-values are 3.84×10⁻⁷ and 4.60×10⁻⁷. Dense regional data are needed because linked variants can all show associations: for example, the conditional H4 leader in the insomnia comparison is rs416223, which carries about 0.726 of conditional SNP probability **given H4**. This is not an unconditional probability of causality and it need not identify the same biological variant under a multi-signal model. The original exact-variant sQTL observations were not included in this expression model.

Brain6's previous trait–trait chromosome 5 coloc showed PP4 0.991 at p12=10⁻⁵ and 0.917 at 10⁻⁶, consistent with a shared signal under its model. The new testis analyses ask a different question: whether a regional expression association shares a signal with each GWAS. Separate pairwise posteriors cannot prove all three traits share one causal variant, establish mediation, or determine direction of regulation in a relevant brain cell. The earlier AlphaGenome AVI and frontal DNase outputs are predictions; they were not used in these fits and do not validate the testis signal.

The original region was selected post hoc and spans one megabase; cis windows may extend beyond it. After filtering, matched variants cover 79.7% and 77.8% of eligible testis QTL records for insomnia and ADHD, respectively; this is not the fraction of all true causal variants covered. The frozen GWAS regional files lack allele frequencies, limiting independent harmonization checks. Source ancestry, overlap, imputation, priors, incomplete variant coverage, and possibly multiple expression or GWAS signals affect inference. The two comparisons share the QTL data and their credible sets are correlated. Multi-signal, ancestry-matched analysis with adequate LD QC, independent cohort replication, brain or developmental molecular data and experimental tests are needed before assigning an effector gene. Neither these results nor the earlier trait–trait result repair Brain6's failed LAVA family gate.

[Protocol](PROTOCOL.md), [posterior table](posteriors.tsv), [conditional SNP table](conditional_h4.tsv), [QTL extraction receipt](qtl_receipt.json), [matching and source-hash receipt](match_receipt.json), [code](extract_qtl.py) and [independent arithmetic check](verification.json) preserve the analysis. All 12 estimated posterior rows were reproduced within 6.9×10⁻¹⁵, which checks arithmetic, not model calibration. GTEx v8 indexed nominal data were read from the [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/); the remote complete files and indices were not fully hashed.

A [common-support sensitivity](common_support/RESULTS.md) finds no variant shared by all three pairwise 95% conditional sets, despite high regional H4 values. The distinct groups have reference r² around 0.70; this does not establish separate causal signals or exclude a shared one.

A subsequent [three-trait HyPrColoc pilot](hyprcoloc_pilot/RESULTS.md), assuming independent study statistics, grouped all three traits under default settings but was sensitive to prior and threshold choices. The result remains exploratory and testis-specific.

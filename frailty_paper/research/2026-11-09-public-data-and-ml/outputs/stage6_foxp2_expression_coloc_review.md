# FOXP2 expression follow-up for the insomnia–frailty locus

## Result

The recent MVP-insomnia/FinnGen-hospital-frailty result supports one possible shared **disease association signal** in the selected FOXP2-region window. It does not identify an effector gene. I therefore separately compared that signal with two public FOXP2 expression-QTL datasets on their shared variants, using the larger gnomAD NFE and FIN reference matrices. These are pairwise comparisons, not a three-trait colocalization.

| Molecular expression dataset | Disease comparison | Variants | PP.H3 | PP.H4 |
|---|---|---:|---:|---:|
| Pituitary FOXP2 | MVP insomnia | 1,191 | 0.969 | 0.030 |
| Pituitary FOXP2 | FinnGen HFRS | 1,191 | 0.993 | 0.0035 |
| Excitatory-neuron FOXP2 | FinnGen HFRS | 720 | 0.905 | 0.089 |

These are model-conditional posteriors for specific signal pairs with p1=p2=1e-4 and p12=1e-5. H3 corresponds to distinct underlying signals and H4 to a shared signal under the model. The results **favor distinct signals in these particular comparisons**, but missing variants, limited coverage, reference and sample uncertainty mean they do not prove that FOXP2 has no causal role. The excitatory-neuron/MVP comparison could not be reported: on those 720 SNPs the MVP fit returned no credible set under the package's default reporting/purity settings. That is not a negative disease or gene association.

With p12 ranging from 1e-6 to 5e-5, pituitary/MVP PP.H4 ranged 0.0031–0.135; pituitary/HFRS PP.H4 ranged 0.00035–0.0173; brain/HFRS PP.H4 ranged 0.0097–0.329. On the same 1,191 pituitary-matched SNPs, the disease-to-disease signal pair retained PP.H4=0.980 at p12=1e-5. The separate comparisons may involve different variants and signal sets; their posteriors cannot be multiplied into a three-trait probability. They do not validate the original score or establish that brain or pituitary expression mediates either phenotype.

## Insertion sensitivity update

The pituitary credible-set insertion was subsequently matched in both disease GWAS and both gnomAD LD references. Adding it represents all four unique source credible-set variants and preserves the modeled preference for distinct pituitary expression versus disease signals. This remains conditional and does not prove FOXP2 is uninvolved. Details and revised results: `stage6_foxp2_pituitary_indel_coloc_review.md`.

## Coverage and limitations

The pituitary source is eQTL Catalogue QTD000301, FOXP2 expression, with 237 samples in release metadata. The analysis aligned 1,191 of the 1,257 disease/reference SNPs. It includes three of four unique variants in the published molecular credible set, covering approximately 0.661 of its total reported PIP; a T>TA insertion with PIP≈0.312 was excluded because the disease/reference setup used SNPs. Missing such a candidate can distort colocalization. The eQTL Catalogue source's per-variant imputation quality was unavailable for this dataset. Its reported genotype/sample ancestry and overlap with each disease study have not been independently established. Z/LD mismatch, effect scaling, small sample size, tissue relevance and out-of-sample LD remain limitations.

The brain source is the public SingleBrain excitatory-neuron full-association release. The source analysis aggregates 983 samples from 757 unique donors across cohorts; an exact gene- and cell-specific effective N is not established here. The 720-SNP intersection came from 1,026 source rows in the extracted window with positive fixed-effect SE and a resolvable effect allele. The reported fixed beta and fixed SE, rather than the portal's random-effect p, were used together. Using working N=757 and N=983 gave brain/HFRS PP.H4 0.0891 and 0.0893 at p12=1e-5; neither proves the correct effective N. Many unmatched or filtered variants remain, and the original rs12536335 is present. No causal claim follows from its marginal FOXP2 association.

Both brain and pituitary analyses use ancestry-oriented gnomAD NFE for the molecular and MVP inputs and FIN for HFRS, not in-sample LD. SuSiE-RSS 0.14.2 used L=10, 95% coverage and residual-variance estimation disabled. All numerical fits converged; Z/LD mismatch estimates were approximately 0.0073, 0.0015 and 0.0018 for pituitary, MVP and HFRS on the 1,191-SNP set, and 0.0018, 0.0010 and 0.0024 for brain, MVP and HFRS on the 720-SNP set. No universal mismatch threshold is established. Convergence does not establish accuracy, and input windows, filtered SNPs, source beta conventions and any sample overlap require scrutiny. The custom gnomAD Hail-index decoder has not yet been checked against Hail itself.

## Reproducibility and source records

Pituitary data and prior single-signal analyses: `stage3_catalogue_foxp2_pituitary_extract_receipt.json`, `stage3_foxp2_pituitary_susie_credible_set.csv`, and `stage5_foxp2_indel_colocalization_review.md`. The SNP-matching and source credible-set coverage are in `stage6_foxp2_pituitary_gnomad_harmonization.json` and `stage6_foxp2_pituitary_gnomad_aligned.csv`; scripts are `stage6_foxp2_pituitary_gnomad_harmonize.py` and `stage6_foxp2_pituitary_gnomad_coloc.R`.

For brain, the targeted source extraction, including BGZF-range provenance, is `stage6_singlebrain_ext_foxp2_regional_receipt.json` and `stage6_singlebrain_ext_foxp2_full_regional.tsv.gz`; scripts `stage6_singlebrain_foxp2_regional_extract.py Ext`, `stage6_foxp2_singlebrain_gnomad_harmonize.py`, `stage6_foxp2_singlebrain_gnomad_coloc.R`, and the N=983 sensitivity script reproduce the analysis. The source's full-file checksum was not verified because only a range was retrieved. Complete result tables and diagnostics follow the `stage6_foxp2_pituitary_gnomad_*` and `stage6_foxp2_singlebrain_gnomad_*` names. Sources: [SingleBrain study](https://www.nature.com/articles/s41588-026-02541-x), [SingleBrain full-statistics release](https://zenodo.org/records/14908182), [eQTL Catalogue data access](https://www.ebi.ac.uk/eqtl/Data_access/), and the disease/reference analysis in `stage6_foxp2_gnomad_matched_coloc_review.md`.

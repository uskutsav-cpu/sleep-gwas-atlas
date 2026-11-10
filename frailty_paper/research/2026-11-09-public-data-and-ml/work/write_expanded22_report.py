import csv,pathlib
p=pathlib.Path('outputs');r=list(csv.DictReader((p/'expanded_22_candidate_risk_gene_table.csv').open()))
s='''# Twenty-two candidate variants: shared-risk follow-up and nearby genes

This is a research shortlist for discussing insomnia–frailty risk, **not a validated risk score**. It contains all 22 of the 109 eligible Song et al. candidates with nominal p<0.05 in both the MVP insomnia and Finnish hospital-frailty (HFRS) GWAS. These markers were selected after prior cross-trait results and follow-up; they are not a prespecified or independent discovery set. All have the same estimated direction in the two follow-up outcomes, but same-sign point estimates do not establish a shared causal variant or mechanism.

| Variant | Counted allele | Nearest gene-body label | FI p | FFS p | Larger of MVP insomnia and HFRS p-values |
|---|---|---|---:|---:|---:|
'''
for x in r:s+=f"| {x['variant']} | {x['counted_allele']} | {x['nearest_gene_body_positional']} | {float(x['FI_p']):.3g} | {float(x['FFS_p']):.3g} | {float(x['joint_max_p']):.3g} |\n"
s+='''
The labels identify the nearest annotated gene body within ±100 kb of the position assumed to be GRCh37. They do not establish an effector gene. Some variants overlap noncoding RNAs or pseudogenes; a neighboring protein-coding gene or a distant regulated gene could be relevant instead. For example, our nearest-body lookup labels rs34290943 as LAMB2, while Song et al. discussed DALRD3 in their broader gene analyses; these represent different mapping approaches, not a confirmed assignment of either gene to this variant.

## What the list can and cannot support

- Six rows meet the numerical max(p) × 124 < 0.05 benchmark in the earlier exploratory screen. The benchmark is not a valid post-selection confirmation without additional assumptions and an independent prespecified analysis.
- Eleven variants have matching point-estimate signs and nominal p<0.05 across the reviewed FI, FFS, UKB insomnia, MVP insomnia and HFRS inputs. These observations are affected by the selection, multiple testing, overlapping or unknown participants and phenotype differences. The FFS effect-allele convention is inferred. Failure of a nominal threshold is not evidence of no effect.
- All 22 specific two-allele contrasts matched biallelic records in the 503-person European 1000 Genomes reference panel. Fourteen rsIDs have more than two alleles in broader Ensembl GRCh37 records, so this check does not settle original study-specific coding or establish that every source used an equivalent contrast.
- Five same-chromosome pairs within 1 Mb were evaluated in the reference panel. For the tested contrasts, their haplotype r² values range from 0.031 to 0.160; the highest was rs1909650–rs9536347. These values are specific to this panel, and low pairwise r² does not establish conditional independence or shared causality. Long-range LD, other variants, and population differences were not exhaustively evaluated.
- The originally prepared 15-marker FI, FFS and insomnia weighted sums are separate technical prototypes, not absolute risk probabilities or evaluated predictors. This 22-marker set was selected using MVP and HFRS results and has **not** been fitted, calibrated, or evaluated against independent outcomes. Do not add FI, FFS, HFRS and insomnia coefficients, which use different scales.

## Recommended next analysis

Freeze a study question and candidate-selection rule before looking at a new validation cohort. Verify alleles and ancestry-matched LD, consider locus-level fine-mapping and colocalization, and link variants to genes with regulatory evidence. For a 15–30-variant score, choose one outcome and time horizon, select and estimate weights in a development set, then evaluate it in independent genotyped participants against clinical baseline, FI and FFS scores. Assess calibration and incremental prediction; perform a sleep/fatigue-excluded FI sensitivity if data permit. No such outcome-based evaluation was possible with the reference genotypes used here.

## Files and provenance

`expanded_22_candidate_risk_gene_table.csv` retains coordinates, the counted/other alleles, signed source effects, p-values, gene annotations, and reference flags. Its SHA-256 is `d500617af1a76f8d7408c26a2b4ba277a3f736335aab7efd3f364459872fc874`. `song124_22_nominal_positional_genes.csv` preserves Ensembl annotation URLs and response hashes. `song124_22_nominal_reference_ld.json` records reference-panel matching and the five local pairwise LD estimates. The input screen and assumptions are documented in `song124_mvp_hfrs_followup.md`; the presumed GRCh37 build of the Song table remains not independently authenticated. No individual-level genotypes or scores were saved.
'''
(p/'expanded_22_candidate_risk_gene_review.md').write_text(s)

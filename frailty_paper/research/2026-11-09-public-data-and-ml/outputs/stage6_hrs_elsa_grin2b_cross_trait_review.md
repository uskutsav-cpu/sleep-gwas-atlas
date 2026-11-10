# HRS–ELSA frailty locus in public insomnia and hospital-frailty results

## Findings

[Mekli et al. (2018)](https://doi.org/10.1371/journal.pone.0207824) reported a frailty-index association for rs7134291 near GRIN2B in HRS and a nominal same-direction association in ELSA. Their published supplement lists the A allele, with HRS beta −0.01021 (p=1.853e−6; 8,227 genotyped individuals in that row) and ELSA beta −0.005013 (p=0.03433; 5,248 individuals). Both outcomes are square-root-transformed FIs; ELSA replication is nominal and the discovery result was selected. This is an independent-cohort *reported frailty locus*, not validation of the project's 15- or 22-marker score.

I checked this exact A/G variant (GRCh38 chr12:13975388) in two public GWAS sources. Both had one exact allele-pair match, with A as effect allele:

| Public outcome | A-allele estimate | P | Interpretation |
|---|---:|---:|---|
| MVP clinical insomnia | OR 1.0052 (source 95% CI 0.9911–1.0195) | 0.4646 | No clear association in this lookup; the interval allows modest positive or negative effects. |
| FinnGen hospital frailty (HFRS) | beta −0.01961 (SE 0.01347; approximate 95% CI −0.0460 to +0.0068) | 0.1453 | Direction agrees with the published FI estimates, but the HFRS estimate is imprecise and the phenotype and scale differ. |

The candidate was selected from published FI results and the insomnia/HFRS lookups are exploratory. Neither p-value establishes a null effect, and comparing significance across traits does not test whether their effects differ. HRS, ELSA, MVP and FinnGen are distinct named studies, but complete participant and relatedness non-overlap has not been verified. The direction of the effects and the exact A/G string match do not establish a shared causal variant, gene mediation, or sleep–frailty colocalization. GRIN2B is a positional gene annotation. The later [Atkins FI study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8441299/) reports that neither of the earlier study's two highlighted variants was associated with FI in its larger UK Biobank analysis; the source statement is additional context, not a new effect estimate computed here.

## Can this validate our score?

The Mekli supplement contains only 32 selected HRS association rows, not full GWAS statistics. None of their rsIDs exactly matches either the project's 15-marker weighted prototype or its separate 22-marker prototype. That does **not** show those candidates were unmeasured or unassociated in HRS; they may be outside the table's selection threshold. The study's HRS and ELSA genotypes are archived under controlled-access mechanisms, so a publicly viewable paper or data catalogue does not provide unrestricted individual-level genotype-and-outcome data. This source cannot validate the exact score, its calibration, or its incremental predictive value.

## Source checks and reproducibility

`stage6_mekli2018_source_audit.py` extracts the published Supplementary Table S3, checks the 32 selected rsIDs against both prototypes, and records its SHA-256. `stage6_hrs_elsa_lead_external_lookup.py` uses NCBI's GRCh38 rs7134291 mapping, then retrieves position/allele matches from GWAS Catalog's harmonized indexed files. `stage6_hrs_elsa_lead_external_lookup.csv` and its JSON receipt retain the source alleles, index hashes and exact counts. The MVP index checksum matched its published checksum; the HFRS harmonized checksum listing contained no checksum values, so that index was hashed locally but not independently verified against a published hash. Neither complete harmonized file was downloaded or checksum-verified in this lookup. The two source file headers were read separately: the MVP statistic is an odds ratio and the HFRS statistic a beta. The original full HFRS source had been checksum-verified in the earlier candidate audit.

The variant's reference A/G pair was verified from NCBI; the 2018 table specifies its A1 and GRCh37 position; a chain-based lift maps that position to the same GRCh38 coordinate on the positive strand. Residual source-strand/representation uncertainty and participant overlap remain, and the selected-locus comparison is not a substitute for colocalization or a full outcome model. No individual records were accessed.

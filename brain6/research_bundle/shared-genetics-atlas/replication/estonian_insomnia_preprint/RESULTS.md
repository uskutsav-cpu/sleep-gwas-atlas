# Estonian clinical-insomnia follow-up: selected candidate lookup

An unreviewed [clinical sleep-disorder genetics preprint](https://doi.org/10.1101/2025.11.21.25340447) reports discovery meta-analysis across UK Biobank, FinnGen, MVP, Mass General Brigham Biobank and All of Us, followed by an Estonian Biobank analysis. Its supplementary Table S5 provides Estonian results for selected discovery leads, not full regional summary statistics. The authors report 12,962 Estonian insomnia cases.

Two of the 425 previously ranked Brain6 variants match an insomnia row by rsID, GRCh38 position and allele pair. Neither is one of the six original Brain6 regional leads, and neither has a nominal p-value below 0.05:

| Candidate | Region | Estonian clinical-insomnia p | Estonian ALT-aligned beta (SE) |
|---|---|---:|---:|
| rs1093463 C>T | Chromosome 3, rs62264764 region | 0.539 | 0.01149 (0.01871) |
| rs1989903 A>G | Chromosome 7, FOXP2 region | 0.465 | 0.00960 (0.01314) |

These are raw, post-selected sleep-side associations; no joint ADHD test was performed. The source table lists 205,165 participants for both rows. Weak evidence is not proof of no effect, and clinical registry insomnia differs from the original questionnaire measure. The paper describes the Estonian cohort as a follow-up cohort, but overlap with every Brain6 discovery dataset and relatives was not independently audited. No original lead matched this selected table, which cannot establish whether those leads were tested or associated in the full Estonian data.

In a subsequent [post-hoc aligned-effect comparison](DIRECTION_PROTOCOL.md), both Estonian estimates have the same sign as the original Brain6 insomnia estimates, but their Wald 95% intervals cross zero: −0.0252 to 0.0482 for rs1093463 and −0.0162 to 0.0354 for rs1989903. Thus these weak estimates do not reliably confirm the sign. The [comparison table](direction_comparison.tsv), [source hashes](direction_receipt.json) and [interval check](direction_verification.json) document the arithmetic; differing definitions, selection and unknown cohort overlap limit comparison of effect sizes.

The preprint reports 35 clinical-insomnia discovery associations, while the downloaded S5 sheet contains **33 rows labeled insomnia**; the difference is unresolved. This reinforces that the sheet should not be treated as a complete regional resource. Discovery includes sources overlapping other project analyses, and no independent ADHD or shared-causal-variant replication follows.

A [FinnGen ADHD follow-up of these two variants](adhd_followup/RESULTS.md) found exact ADHD records for both. One ADHD p-value was nominal, but neither variant met even a nominal two-trait benchmark because its Estonian insomnia p-value was much larger; this comparison is post-selected and does not establish independent replication.

The [protocol](PROTOCOL.md), [matched rows](matched_candidates.tsv), [source receipt](receipt.json), [lookup](lookup.py) and [independent row verification](verification.json) document the analysis. The supplementary archive and workbook passed ZIP integrity checks; the workbook MD5 matches the Europe PMC full-text XML metadata. The two-candidate extension was made after inspection and is not a prespecified replication family. Failed LAVA gates remain unchanged.

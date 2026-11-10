# FI versus physical-frailty classification in a US veteran cohort

A large [Million Veteran Program study](https://doi.org/10.1161/JAHA.123.033111) measured a questionnaire-based 36-item FI and a modified three-item Study of Osteoporotic Fractures (SOF) physical-frailty measure in the same participants. Reconstructing the authors' Table 1 gives **φ = 0.268** for their binary frail/not-frail classifications among 186,692 participants with both classifications. This is a measured-score classification association, not a genetic correlation, and SOF is **not** the five-component FFS used in the project's Ye GWAS.

The published SOF measure is self-reported and uses weight loss, reduced energy, and a transfer-assistance question as a proxy for a chair-stand test. Its categories are robust (0), pre-frail (1), and frail (2–3). For the FI, the authors define frail as >0.20; I grouped robust and pre-frail together for both measures solely for the binary comparison.

| | Modified SOF frail | Modified SOF not frail | Total |
|---|---:|---:|---:|
| MVP FI frail | 11,923 | 41,897 | 53,820 |
| MVP FI not frail | 6,168 | 126,704 | 132,872 |
| Total | 18,091 | 168,601 | 186,692 |

Binary agreement was **74.3%**, Cohen's κ **0.218**, and the proportion of those identified as frail by either measure who were identified by both was **19.9%**. These are derived from the published table, not new individual-level analyses. The different thresholds and prevalences strongly affect these statistics. Overall agreement includes the large group classified as not frail by both. The counts do not imply one measure is the correct diagnosis, nor compare predictive performance. The parent study includes 190,688 participants, so 3,996 are absent from the SOF cross-tabulation; the reasons and characteristics of those missing participants have not been checked.

For context, the same study's 31-deficit claims-based VA-FI and questionnaire FI binary classifications yield φ=0.409 and κ=0.407 in 190,688 participants. That is a separate comparison using a different sample denominator and different phenotype/thresholds; no paired test was done to show it differs from the SOF association. The cohort was about 94% male and 85% self-identified White, so generalization to other groups needs evaluation. Neither the SOF comparison nor the VA-FI comparison estimates insomnia–frailty association or genetic risk.

## Relevance to score validation

The study shows that multiple frailty measures and longitudinal outcomes have been constructed in a large MVP cohort. It does **not** establish that the same analyzed participants' genotype and sleep data are accessible to this project, that a full FFS is available, or that our score predicts their outcomes. The study used a modified SOF, and its reduced-energy component overlaps conceptually with fatigue. The detailed FI and VA-FI sleep-item and code-level overlap must be checked.

The project's expanded 22-candidate list was selected using an MVP insomnia GWAS, so MVP participants are **not automatically independent** of selection or sleep discovery. Exact participant and relatedness overlap must be checked; for a clean 22-marker external assessment, use a nonoverlapping cohort or redo discovery and selection after excluding the validation individuals and close relatives. The original 15-marker prototype was selected using Song/UKB/FI data rather than MVP follow-up, but its independence from any particular MVP sample also requires verification. Approval and access to genotype-linked frailty and sleep measurements would be required before scoring; no data request or message was sent.

## Reproducibility

`mvp_frailty_agreement.json` retains the two transcribed 3×3 count matrices, collapsed 2×2 counts and formulas' resulting statistics. Counts were checked against published Table 1 and the [authors' accessible preprint](https://www.medrxiv.org/content/10.1101/2023.10.18.23297231v1.full.pdf); the published article describes the 36-item FI. The estimates are derived from aggregate cross-tabulation, with no individual-level data or standard errors computed.

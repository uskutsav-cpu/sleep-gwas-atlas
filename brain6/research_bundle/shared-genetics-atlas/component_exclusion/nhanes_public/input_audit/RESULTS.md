# NHANES input readiness and walking-question skip pattern

The public 2011–2012 and 2013–2014 demographic, depression-screen, physical-function and grip files were retrieved from the CDC and parsed. All selected files had unique participant IDs and contained the variables needed for this **initial** check. This is not a completed 49-item FI, adapted Fried score, phenotype correlation or genetic analysis.

| Unweighted count among participants aged 60+ | 2011–2012 | 2013–2014 |
| --- | ---: | ---: |
| Demographic records | 1,791 | 1,841 |
| Valid sleep and fatigue responses and valid quarter-mile and lifting responses | 1,221 | 1,293 |
| Those four responses plus at least one positive grip reading | 1,139 | 1,217 |
| Reported difficulty walking without equipment | 383 | 411 |
| Of those, quarter-mile question missing | 383 | 411 |

The [CDC codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/PFQ_G.htm) states that the quarter-mile question is not displayed when the participant reports difficulty walking without equipment. The corresponding [2013–2014 codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/PFQ_H.htm) should also be used when completing cycle-specific coding. A complete-case analysis of that question would exclude many participants with mobility limitations; treating their missing answers as no difficulty would be wrong. How to score or impute the skipped item needs an explicit, source-justified rule and sensitivity analysis. The published adapted Fried study's methods alone do not establish how it handled all skipped responses.

Counts are file-level sample counts, **not population estimates**; no survey weights or variance estimates were applied. The four-response subset does not establish complete FI or Fried-score coverage. Some additional excluded responses may reflect refusal, nonresponse or `do not do this activity`; the valid-code definitions in the [protocol](PROTOCOL.md) and [audit](audit.R) are intentionally narrow and are not a finalized phenotype recipe. Grip availability means at least one positive trial, not a validated weakness score. The two instruments reuse physical questions, and their definitions differ from the original project’s UK Biobank GWAS measures.

The eight original public XPT files and [download receipt](download_receipt.json) are on the SSD-backed `/Volumes/Brain6Sim/nhanes-public-2011-14-20261028`; the receipt records CDC URLs, file sizes and SHA-256. The XPT files passed a header check and `foreign::read.xport` parsing, but no independent publisher checksum was available in this check. [File counts](files.tsv) and [all readiness counts](counts.tsv) are saved here. A separate [match-based check](verify.R) reproduced the sleep, fatigue, physical-question and skip counts; it did not independently check the grip or survey-design rows. No participant-level derived data were saved to the atlas. Further documentation, the remaining FI and Fried components, matching rules and a survey-design analysis are needed before reporting correlations.

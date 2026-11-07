# Clinical ADHD, insomnia and hospital frailty at selected markers

This post-hoc join uses the seven markers in the existing FinnGen clinical ADHD lookup, the saved MVP insomnia results and the existing Finnish hospital-frailty extraction. A subsequent complete-source extraction obtained the missing hospital-frailty row for rs2582894. Six markers have all three p-values; one remains not estimable. **None of the six meets the descriptive factor-of-seven max-p benchmark of 0.05.**

| Marker | MVP insomnia p | FinnGen clinical ADHD p | FinnGen hospital-frailty p | 7 × largest p, capped at 1 |
|---|---:|---:|---:|---:|
| rs2582894 | 0.000933 | 0.0653 | 0.0000509 | 0.457 |
| rs7160630 | 0.000791 | 0.0133 | 0.0775 | 0.543 |
| rs62264767 | 0.0000499 | 0.0261 | 4.63×10⁻⁷ | 0.183 |
| rs2431108 | 0.0895 | 0.000108 | 0.00110 | 0.627 |
| rs6452785 | 0.0000310 | No indexed exact row | No matching rsID in HFRS scan | Not estimable |
| rs1476535 | 1.20×10⁻⁷ | 0.367 | 1.46×10⁻⁸ | 1.000 |
| rs62264764 | 0.0000492 | 0.0279 | 4.52×10⁻⁷ | 0.195 |

rs62264764 and rs62264767 are highly linked and the first was added separately after previous examination; they are not independent loci. The maximum p-value is a conservative intersection–union p-value for marginal association with all three outcomes if each source p-value is valid. It does not require the tests to be independent, but it cannot establish a common causal variant, gene, mediation or interaction. Multiplication by seven is a numerical benchmark for this selected set and does not correct the earlier, broader project search. Non-crossing is not proof of no shared biology.

The clinical ADHD phenotype differs from the quantitative childhood symptom measure in the earlier [425-candidate screen](../../variant_context/brain6_three_trait_candidate/RESULTS.md). MVP and FinnGen are separate named sources, but exact overlap with discovery cohorts and relatedness are unverified; both clinical ADHD and hospital frailty here come from FinnGen. HFRS is hospital-ascertained and contains a somnolence-related diagnostic category, not a sleep/exhaustion-excluded physical-frailty composite. Effects and uncertainty differ between sources, and MVP effect directions derive from rounded odds ratios. Clinical ADHD sample size, Finnish LD, quality and the incomplete indexed source checks remain as described in the [clinical ADHD source audit](../finngen_adhd_raw_lookup/RESULTS.md). The missing HFRS rsID and ADHD indexed record for rs6452785 are not treated as null associations.

The [protocol](PROTOCOL.md), [initial five-match results](results.tsv), [initial join code](run.py) and [input hash receipt](receipt.json) retain the original calculation. The earlier MVP/HFRS source files were read fully and hash checked, while the clinical source was checked by indexed extraction and a separate whole-source compressed-byte hash, without complete decompression. These checks do not establish statistical calibration or causal sharing.

The [subsequent full-source scan](complete_missing.py) read all 21,294,560 HFRS rows, reproduced the previously recorded MD5 and found the unique coordinate- and allele-matched rs2582894 result. No rs6452785 rsID was found in that scan; this is not evidence of no association. The [completed table](completed_results.tsv), [matched row](rs2582894_completion.tsv), and [receipt](completion_receipt.json) preserve the addition. The [table builder](build_completed.py) reproduces the six-match table, verifies its input hashes and checks every max-p calculation; its [receipt](completed_table_receipt.json) records the output hash. The original seven-slot benchmark and conclusions remain unchanged.

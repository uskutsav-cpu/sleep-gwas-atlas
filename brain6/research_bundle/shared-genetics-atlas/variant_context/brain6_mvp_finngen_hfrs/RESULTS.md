# Selected Brain6 variants in MVP insomnia and Finnish hospital frailty

A later [three-phenotype follow-up](../brain6_three_trait_candidate/RESULTS.md) checks the 425 candidates against childhood ADHD symptoms. It reports a nominal rs4916723 symptom result and no candidate crossing its full three-phenotype benchmark.

Among the 425 previously ranked Brain6 variants, **419** matched both the MVP European-ancestry insomnia GWAS and FinnGen hospital-frailty GWAS by rsID, GRCh38 position and non-palindromic allele pair. **155 variants** have `max(p_insomnia, p_HFRS) × 425 < 0.05`: 92 in the selected chromosome 3 region, 36 on chromosome 11 and 27 on chromosome 7. All 155 have the same estimated direction for the selected ALT allele across the two outcomes. They are correlated SNPs in three selected regions, **not 155 independent loci**. Their direction agreement is descriptive and was not itself tested.

| Brain6 candidate | MVP insomnia p | Finnish HFRS p | max-p × 425 | ALT-aligned directions | Childhood ADHD-symptom p |
| --- | ---: | ---: | ---: | --- | ---: |
| rs62264764, chr3, G>A | 4.918×10⁻⁵ | 4.522×10⁻⁷ | 0.0209 | Both negative for A | 0.437 |
| rs62264767, chr3 | 4.993×10⁻⁵ | 4.626×10⁻⁷ | 0.0212 | Both negative for ALT | 0.478 |
| rs2431108, chr5 | 0.0895 | 0.00110 | 1.000 | Both positive for ALT | 0.00600 |
| rs77960, chr5 | 0.0860 | 0.000967 | 1.000 | Both positive for ALT | 0.00606 |
| rs4916723, chr5 | 3.345×10⁻⁵ | 0.00919 | 1.000 | Both negative for ALT | Not evaluated in the cited follow-up |
| rs1476535, chr7, C>T | 1.201×10⁻⁷ | 1.456×10⁻⁸ | 0.0000510 | Both positive for T | 0.660 |

The ADHD-symptom values are previously obtained separate follow-up results, not part of the two-trait max-p test. rs62264767, rs2431108 and rs1476535 are from the [six-lead follow-up](../../replication/six_lead_joint/RESULTS.md); rs62264764 and rs77960 are from the [three-variant symptom lookup](../../replication/childhood_symptoms/RESULTS.md). Their quantitative childhood symptoms are not diagnostic ADHD, and weak p-values do not prove no effect. In particular, the rs62264764 and rs1476535 sleep–HFRS support does **not** establish independent replication across sleep, frailty and ADHD.

The max-p is an intersection–union test for marginal association with both phenotypes if the marginal p-values are valid, without requiring the two source tests to be independent. The 425-slot factor includes all ranked candidates, even those unavailable in a source. However, this family and its regions were selected after prior analyses, and some source findings and proxies were previously inspected; the displayed benchmark is not a confirmatory correction for the broader search. It does not test one shared causal variant, a causal gene or a mediation pathway. In particular the chr11 region has low H4 in the earlier single-signal insomnia–ADHD model; the new marginal sleep–HFRS results do not override that evidence. LD and multiple causal signals require appropriate fine-mapping and colocalization, with independently characterized cohorts and reference ancestry.

The direct rs62264764 results add sleep–frailty context to the earlier [UKB HFRS screen](../brain6_ukb_hfrs/RESULTS.md), where its original-HFRS p=0.000127 did not meet that screen's different 850-slot benchmark. The chr7 lead is also associated in that UKB HFRS screen. The candidate and source selection, phenotype differences, exact participant/relative overlap, and regional LD remain limitations. No result validates the AlphaGenome prediction, establishes LSAMP, LINC03051 or FOXP2 as a causal target, or resolves the failed LAVA gate.

## Provenance and checks

We read the full 21,294,560-row Finnish HFRS file and 19,703,815-row MVP insomnia gzip, matched their earlier recorded MD5 and SHA-256 hashes respectively, and read the gzip through completion. HFRS yielded 420 exact candidate matches and MVP 424; six candidates were unavailable in at least one source and remain in the 425-slot denominator. Allele signs were aligned to the candidate ALT; MVP uses the log of a rounded source odds ratio, so it should not be treated as a full-precision effect. The source-p values, not p-values reconstructed from those odds ratios, drive the test. Source positions and alleles are consistent with the candidate crosswalk, but original GWAS QC and source model validity are not established by matching.

The complete result table is `/Volumes/Extreme SSD/brain6-mvp-finngen-hfrs-results-20261028.tsv`, the reproducible script is `/Volumes/Extreme SSD/brain6-mvp-finngen-hfrs-20261028.py`, and source and output hashes are in `/Volumes/Extreme SSD/brain6-mvp-finngen-hfrs-receipt-20261028.json`. The [protocol](PROTOCOL.md) records the test before this 425-variant extraction. A separate arithmetic check reproduced all 419 paired max-p values, 425-factor adjustments and direction flags.

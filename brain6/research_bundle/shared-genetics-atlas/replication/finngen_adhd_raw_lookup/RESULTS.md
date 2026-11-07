# Public clinical ADHD follow-up of Brain6 leads

The public **FinnGen R13 F5_ADHD** indexed release contains exact coordinate, rsID and allele-pair matches for five of the six original Brain6 leads and for the separately selected rs62264764. The chromosome 5 lead **rs2431108** has C-allele beta +0.0825 (SE 0.0213), p = 0.000108. Its effect has the same direction as the frozen discovery ADHD input (+0.0730 on C). This is encouraging clinical ADHD-side follow-up at a previously selected locus, not a prospective or selection-adjusted replication claim.

| Variant | GRCh38 coordinate | Aligned allele | ADHD beta (SE) | ADHD p | Direction vs discovery ADHD |
|---|---:|---|---:|---:|---|
| rs2582894 | 11:28554714 | G | −0.0607 (0.0329) | 0.0653 | Same |
| rs7160630 | 14:98171142 | G | +0.0687 (0.0277) | 0.0133 | Same |
| rs62264767 | 3:117923158 | C | −0.0530 (0.0238) | 0.0261 | Same |
| rs2431108 | 5:104612267 | C | +0.0825 (0.0213) | 0.000108 | Same |
| rs6452785 | 5:88389683 | T | Not available | Not available | Not estimable |
| rs1476535 | 7:114430980 | C | −0.0177 (0.0196) | 0.367 | Same, imprecise |
| rs62264764 | 3:117920728 | A | −0.0524 (0.0238) | 0.0279 | Separate selected marker; not checked here |

The index returned no row at the rs6452785 coordinate; absence is missing evidence, not a null association. The separately selected rs62264764 is in very strong LD with rs62264767 in the previously examined European reference, so these are not independent regional confirmations. Its nominal ADHD association does not establish a causal variant or link it to LSAMP. [FinnGen's format documentation](https://finngen.gitbook.io/documentation/r10/data-description) defines beta relative to the alternate allele and on the log-odds scale for these case–control results. Effect and SE are on that scale; no absolute risk or clinical application follows.

For rs2431108, the [prior FinnGen R13 insomnia lookup](../six_lead_joint/RESULTS.md) reports C-allele beta +0.0281 and p = 0.000124, also matching its discovery sleep direction. The larger of the two follow-up p-values is 0.000124; multiplied by six it is approximately 0.000744 as a **descriptive six-lead benchmark**. This calculation is not a confirmatory adjusted replication: this candidate and the insomnia value were already inspected, the original discovery selection was broader than six, and both follow-up traits come from FinnGen with unverified participant/relative overlap. Max-p is an intersection–union p-value under valid marginal tests without requiring independence, but it cannot establish a common causal allele or correct the retrospective selection. Registry insomnia differs from the discovery questionnaire phenotype.

The historical [source audit](../finngen_adhd_archived_munged/historical_source_audit.md) reports 5,559 cases and 489,493 controls and describes FinnGen as independent of the discovery ADHD cohorts by published cohort descriptions; this is not a participant-level independence or relatedness audit, and the Finnish sample, ascertainment and LD differ. The rs1476535 raw-source Z (−0.0176867/0.0196199 ≈ −0.9015 on C) is consistent with the previous archived processed row, but a single match does not prove its provenance.

## Reproducibility and limits

[Protocol](PROTOCOL.md), [indexed query](indexed_scan.py), [all records](indexed_results.tsv), [receipt](indexed_receipt.json), [publisher index](../../PUBLICATION.md#external-dependencies) (external target: `source.tbi`) and [source response metadata](response_headers.txt) are retained. The index's full MD5 matches its publisher ETag. The indexed query fetched small source ranges. A separate [whole-source hash pass](whole_source_hash_receipt.json) streamed all 802,195,177 compressed bytes without retaining them: the MD5 matched the historical lock, and the calculated SHA-256 is recorded. The pass did not decompress and validate every gzip member or prove the older munged file's processing lineage. The initial full-source parser failed to handle concatenated BGZF members and was stopped; no result from it is used. The frozen regional input files' hashes and discovery directions were checked separately. These findings do not repair the failed LAVA quality gates, validate AlphaGenome predictions or establish gene mediation.

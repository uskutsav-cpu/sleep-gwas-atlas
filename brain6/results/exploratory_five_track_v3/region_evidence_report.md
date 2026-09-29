# Five-track PLACO exploratory 20-region evidence — current diagnostic ledger

**EXPLORATORY. Zero LAVA-confirmed regions and zero final shared-locus tiers.**

All 25 pair-specific candidates remain in the 20 coordinate-union
groups below. Across 2,686 candidate rows, 2,656 have valid signed
UKB-reference LD/effect orientation, six retain raw |r| > 1 as QC
exceptions, and 24 lack an exact frozen-reference match. All 27 lead
rows pass signed LD. Direction agreement is a post-selection
descriptive check and cannot establish a shared causal effect.

| Geographic region | Pairs | Leads | Candidate rows | Valid signed LD | Invalid raw r | Pair-level global replication | Positional genes (four-pair / B) |
|---|---|---|---:|---:|---:|---|---|
| chr2:57487593-58487593 | longsleep__scz | rs11682175 | 41 | 41 | 0 | NONE | CTD-2026C7.1 / NOT_APPLICABLE |
| chr2:135040546-136040546 | longsleep__parkinson | rs7599054 | 2 | 1 | 1 | NONE | CCNT2-AS1 / NOT_APPLICABLE |
| chr3:52219816-53219816 | longsleep__scz | rs10933 | 3 | 3 | 0 | NONE | PBRM1 / NOT_APPLICABLE |
| chr3:117142005-118142005 | insomnia__adhd | rs62264767 | 92 | 90 | 2 | insomnia__adhd | NOT_ANNOTATED / LSAMP |
| chr4:169728549-170728549 | longsleep__scz | rs4692709 | 1 | 1 | 0 | NONE | SH3RF1 / NOT_APPLICABLE |
| chr5:87185500-88185500 | insomnia__adhd | rs6452785 | 43 | 43 | 0 | insomnia__adhd | NOT_ANNOTATED / TMEM161B-AS1 |
| chr5:103447968-104481726 | insomnia__adhd;insomnia__mdd | rs2431108;rs40465 | 177 | 177 | 0 | insomnia__adhd | RP11-6N13.1 / RP11-6N13.1 |
| chr6:104869510-105869510 | insomnia__mdd | rs4946651 | 5 | 5 | 0 | NONE | LINC00577 / NOT_APPLICABLE |
| chr7:104429267-105429267 | longsleep__scz | rs6466056 | 198 | 198 | 0 | NONE | SRPK2 / NOT_APPLICABLE |
| chr7:113559156-114571035 | insomnia__adhd;insomnia__mdd | rs1476535;rs2894699 | 78 | 78 | 0 | insomnia__adhd | AC073626.2 / FOXP2 |
| chr11:28076261-29076261 | insomnia__adhd | rs2582894 | 1 | 1 | 0 | insomnia__adhd | NOT_ANNOTATED / RP11-22P4.1 |
| chr11:61081764-62081764 | longsleep__bipolar | rs174560 | 19 | 19 | 0 | NONE | FADS1 / NOT_APPLICABLE |
| chr11:112877488-113877488 | longsleep__scz | rs7111031 | 28 | 27 | 1 | NONE | DRD2 / NOT_APPLICABLE |
| chr13:53365563-54365563 | insomnia__mdd | rs9536383 | 49 | 49 | 0 | NONE | RN7SL618P / NOT_APPLICABLE |
| chr14:41646572-42646572 | insomnia__mdd | rs2121708 | 12 | 12 | 0 | NONE | LRFN5 / NOT_APPLICABLE |
| chr14:98137479-99137479 | insomnia__adhd | rs7160630 | 36 | 36 | 0 | insomnia__adhd | NOT_ANNOTATED / RP11-61O1.1 |
| chr15:84700520-85700520 | longsleep__bipolar;longsleep__scz | rs1051168 | 11 | 11 | 0 | NONE | NMB / NOT_APPLICABLE |
| chr17:1335482-2335482 | longsleep__bipolar | rs4790841 | 1 | 1 | 0 | NONE | RTN4RL1 / NOT_APPLICABLE |
| chr17:43049526-45343136 | longsleep__parkinson;longsleep__scz | rs17631676;rs199443;rs199528;rs62054437;rs62073157 | 1845 | 1843 | 2 | NONE | NSF;PLEKHM1;RP11-105N13.4;RP11-259G18.3;WNT3 / NOT_APPLICABLE |
| chr18:52604019-53604019 | insomnia__mdd | rs4468713 | 20 | 20 | 0 | NONE | TCF4 / NOT_APPLICABLE |

The detailed table adds exact 1000G genotype coverage, five measured
cross-pair overlap edges, source-specific GWAS Catalog context, and
explicit empty evidence states. Those five independent-genotype LD
comparisons have r² from 0.755 to 1.000; this corroborates LD among
overlapping leads without establishing independent signals.
The 24 reference-unmatched rows are preserved in the v1 unassigned
table, outside the 20 coordinate groups.

Only insomnia-ADHD has archived **global pair-level** directional
replication in FinnGen; no region has locus-specific replication.
Neither a GWAS Catalog exact-rsID record nor a nearby gene is
independent biological support. No region has validated SuSiE,
trait-trait coloc, eQTL/sQTL coloc, regulatory, tissue/cell-type, or
pathway support. Thus no biologically supported region can be ranked
as a validated shared locus. The multi-pair groups on chromosomes
5, 7, 15, and 17 remain descriptive hypotheses, without a formal
order or evidence-tier promotion.

Canonical LAVA v3 remains `FAILED_QC_NOT_PROMOTED`: 3,720/17,465
NOT_RUN versus a frozen maximum of 873. Local-rg claims, final
region tiers, and confirmatory cross-layer claims are `BLOCKED_LAVA`.
Any further molecular result in this branch stays exploratory and
cannot enter the frozen tier decision.

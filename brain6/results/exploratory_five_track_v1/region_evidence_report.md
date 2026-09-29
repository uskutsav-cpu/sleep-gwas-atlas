# Five-track PLACO exploratory region evidence

**EXPLORATORY; no LAVA-confirmed region or final shared-locus tier.**

The five-track PLACO family has 25 pair-specific candidate intervals. Their
coordinate union comprises 20 geographic groups; neither count is an
independent causal-signal count. Of 2,686 candidate variant rows, 2,662
are assigned to these intervals and 24 lack an exact frozen-reference
match and remain in `unassigned_reference_unmatched_candidates.tsv`.

| Geographic region | Pairs | Leads | Candidate rows | Exact 1000G genotypes | Prior catalog records | Positional gene context |
|---|---|---|---:|---:|---:|---|
| chr2:57487593-58487593 | longsleep__scz | rs11682175 | 41 | 41 | 10 | CTD-2026C7.1 |
| chr2:135040546-136040546 | longsleep__parkinson | rs7599054 | 2 | 2 | 0 | CCNT2-AS1 |
| chr3:52219816-53219816 | longsleep__scz | rs10933 | 3 | 3 | 0 | PBRM1 |
| chr3:117142005-118142005 | insomnia__adhd | rs62264767 | 92 | 92 | 0 | NOT_ANNOTATED |
| chr4:169728549-170728549 | longsleep__scz | rs4692709 | 1 | 1 | 0 | SH3RF1 |
| chr5:87185500-88185500 | insomnia__adhd | rs6452785 | 43 | 43 | 0 | NOT_ANNOTATED |
| chr5:103447968-104481726 | insomnia__adhd;insomnia__mdd | rs2431108;rs40465 | 177 | 173 | 5 | RP11-6N13.1 |
| chr6:104869510-105869510 | insomnia__mdd | rs4946651 | 5 | 5 | 0 | LINC00577 |
| chr7:104429267-105429267 | longsleep__scz | rs6466056 | 198 | 198 | 0 | SRPK2 |
| chr7:113559156-114571035 | insomnia__adhd;insomnia__mdd | rs1476535;rs2894699 | 78 | 77 | 2 | AC073626.2 |
| chr11:28076261-29076261 | insomnia__adhd | rs2582894 | 1 | 1 | 0 | NOT_ANNOTATED |
| chr11:61081764-62081764 | longsleep__bipolar | rs174560 | 19 | 19 | 2 | FADS1 |
| chr11:112877488-113877488 | longsleep__scz | rs7111031 | 28 | 28 | 0 | DRD2 |
| chr13:53365563-54365563 | insomnia__mdd | rs9536383 | 49 | 49 | 3 | RN7SL618P |
| chr14:41646572-42646572 | insomnia__mdd | rs2121708 | 12 | 12 | 0 | LRFN5 |
| chr14:98137479-99137479 | insomnia__adhd | rs7160630 | 36 | 36 | 0 | NOT_ANNOTATED |
| chr15:84700520-85700520 | longsleep__bipolar;longsleep__scz | rs1051168 | 11 | 11 | 0 | NMB |
| chr17:1335482-2335482 | longsleep__bipolar | rs4790841 | 1 | 1 | 1 | RTN4RL1 |
| chr17:43049526-45343136 | longsleep__parkinson;longsleep__scz | rs17631676;rs199443;rs199528;rs62054437;rs62073157 | 1845 | 1838 | 1 | NSF;PLEKHM1;RP11-105N13.4;RP11-259G18.3;WNT3 |
| chr18:52604019-53604019 | insomnia__mdd | rs4468713 | 20 | 20 | 0 | TCF4 |

`geographic_region_evidence.tsv` records the full pair membership,
candidate count, diagnostic genotype/LD coverage, signed-LD reuse status,
pair-level global replication, prior GWAS Catalog study IDs, and the
method-specific empty evidence states for every region. The six
insomnia-ADHD candidate regions inherit a **pair-level** FinnGen global
directional replication label; no region has locus-specific replication.
The existing signed-LD diagnostic applies only where its exact
four-pair candidate-to-lead assignments match this five-track set.
It is not reusable for three loci with changed clumps or any of the
six admitted B loci. Independent 1000G LD is diagnostic and does not
replace the frozen UKB locus definition.

No region has validated SuSiE, trait-trait coloc, eQTL/sQTL coloc,
regulatory, tissue/cell-type, or pathway support in this ledger.
The nearest-gene names are descriptive positions. Curated catalog
matches may reuse discovery cohorts and do not establish independent
replication. Accordingly, there is **no biologically supported region
that can be ranked as a validated shared locus**. The four multi-pair
groups on chromosomes 5, 7, 15, and 17 are descriptive overlap
hypotheses, retained with all of their pair-specific members.

The canonical LAVA v3 family remains `FAILED_QC_NOT_PROMOTED`
(3,720/17,465 NOT_RUN, frozen maximum 873). All local-rg claims,
confirmatory region tiers, and dependent cross-layer categories are
`BLOCKED_LAVA`. Biological annotations in this exploratory branch
cannot enter the frozen tier rule.

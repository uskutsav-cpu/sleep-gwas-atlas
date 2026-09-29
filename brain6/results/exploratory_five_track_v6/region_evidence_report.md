# Brain6 five-track exploratory evidence across 20 geographic regions

**EXPLORATORY. Canonical LAVA v3 remains FAILED_QC_NOT_PROMOTED; no final shared-locus tier is assigned.**

The source-bound table retains all 25 pair-specific PLACO candidate loci, 20 coordinate
groups, 2,686 candidate rows, 27 lead rows, and all prior evidence columns. Pair-level
FinnGen directional replication applies only to insomnia–ADHD globally; it is not
locus-specific replication for any of its six rows below. The GTEx columns distinguish
nominal lead eQTL associations from exact PLACO-candidate overlap with published
molecular credible sets. Neither is GWAS fine-mapping or colocalization.

| Geographic region | Pair(s) | Lead(s) | Valid signed LD / candidates | Global pair replication | Exact lead eQTL queries with rows | Candidates in molecular QTL sets | Highest molecular eQTL / sQTL PIP | Exact lead regulatory features | Prior trait-relevant Catalog rows |
|---|---|---|---:|---|---:|---:|---|---:|---:|
| chr2:57487593-58487593 | longsleep__scz | rs11682175 | 41/41 | NONE | 2/2 | 0 | none / none | 1 | 10 |
| chr2:135040546-136040546 | longsleep__parkinson | rs7599054 | 1/2 | NONE | 2/2 | 0 | none / none | 2 | 0 |
| chr3:52219816-53219816 | longsleep__scz | rs10933 | 3/3 | NONE | 2/2 | 0 | none / none | 1 | 0 |
| chr3:117142005-118142005 | insomnia__adhd | rs62264767 | 90/92 | insomnia__adhd | 2/2 | 0 | none / none | 0 | 0 |
| chr4:169728549-170728549 | longsleep__scz | rs4692709 | 1/1 | NONE | 2/2 | 0 | none / none | 0 | 0 |
| chr5:87185500-88185500 | insomnia__adhd | rs6452785 | 43/43 | insomnia__adhd | 0/2 | 0 | none / none | 0 | 0 |
| chr5:103447968-104481726 | insomnia__adhd;insomnia__mdd | rs2431108;rs40465 | 177/177 | insomnia__adhd | 4/4 | 0 | none / none | 0 | 5 |
| chr6:104869510-105869510 | insomnia__mdd | rs4946651 | 5/5 | NONE | 2/2 | 0 | none / none | 0 | 0 |
| chr7:104429267-105429267 | longsleep__scz | rs6466056 | 198/198 | NONE | 2/2 | 0 | none / none | 0 | 0 |
| chr7:113559156-114571035 | insomnia__adhd;insomnia__mdd | rs1476535;rs2894699 | 78/78 | insomnia__adhd | 4/4 | 0 | none / none | 0 | 2 |
| chr11:28076261-29076261 | insomnia__adhd | rs2582894 | 1/1 | insomnia__adhd | 2/2 | 0 | none / none | 0 | 0 |
| chr11:61081764-62081764 | longsleep__bipolar | rs174560 | 19/19 | NONE | 2/2 | 5 | 0.199166133140605 / none | 1 | 2 |
| chr11:112877488-113877488 | longsleep__scz | rs7111031 | 27/28 | NONE | 2/2 | 0 | none / none | 0 | 0 |
| chr13:53365563-54365563 | insomnia__mdd | rs9536383 | 49/49 | NONE | 2/2 | 0 | none / none | 0 | 3 |
| chr14:41646572-42646572 | insomnia__mdd | rs2121708 | 12/12 | NONE | 0/2 | 0 | none / none | 0 | 0 |
| chr14:98137479-99137479 | insomnia__adhd | rs7160630 | 36/36 | insomnia__adhd | 2/2 | 0 | none / none | 0 | 0 |
| chr15:84700520-85700520 | longsleep__bipolar;longsleep__scz | rs1051168 | 11/11 | NONE | 2/2 | 3 | 0.00220713919388604 / none | 0 | 0 |
| chr17:1335482-2335482 | longsleep__bipolar | rs4790841 | 1/1 | NONE | 0/2 | 0 | none / none | 1 | 1 |
| chr17:43049526-45343136 | longsleep__parkinson;longsleep__scz | rs17631676;rs199443;rs199528;rs62054437;rs62073157 | 1843/1845 | NONE | 8/10 | 1686 | 0.431615618148346 / 0.0332802036973561 | 0 | 1 |
| chr18:52604019-53604019 | insomnia__mdd | rs4468713 | 20/20 | NONE | 2/2 | 0 | none / none | 0 | 0 |

The separate TSV carries positional genes, source-specific Catalog studies, five
independent-genotype cross-pair LD edges, all signed-LD QC exceptions, and every
method state. Six candidate rows have out-of-range raw signed correlations and 24
have no exact frozen-reference allele; all 27 lead rows pass signed LD. In the
regulatory lookup, six feature rows overlap five exact leads; these are coordinate
annotations only. An empty exact-lead query cannot exclude a regulatory feature
or QTL elsewhere in the region.

## Current interpretation

The chr11:61.1–62.1 Mb long-sleep–bipolar region has five candidate variants in
one published brain eQTL credible set, a molecular maximum PIP of 0.199, and an
enhancer coordinate overlap at its lead. The chr17:43.0–45.3 Mb long-sleep–PD/SCZ
group has extensive candidate overlap with published brain expression and splicing
credible sets; it also contains 1,845 PLACO candidate rows, so raw overlap count
is particularly sensitive to candidate density. These are descriptive molecular
hypotheses, not ranked or validated shared loci. The chr15 group also has candidate
eQTL-set overlap, with low molecular PIP in this panel. Nearby gene labels and
Catalog rows remain positional/prior-study context; they do not establish causality.

**No region has locus-specific independent replication, validated GWAS SuSiE
fine-mapping, trait–trait coloc, trait–eQTL/sQTL coloc, tissue/cell-type
enrichment, or pathway support.** The available nominal eQTL and molecular
credible-set observations do not fill those missing method states. LAVA local-rg,
final evidence tiers, and confirmatory cross-layer claims remain BLOCKED_LAVA.

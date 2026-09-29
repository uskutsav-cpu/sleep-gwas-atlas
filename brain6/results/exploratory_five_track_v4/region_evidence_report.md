# Five-track PLACO exploratory 20-region evidence — lead eQTL context added

**EXPLORATORY. Zero LAVA-confirmed regions and zero final shared-locus tiers.**

All 25 pair-specific candidates, 27 lead rows, and 20 coordinate-union groups are retained.
Exact lead GRCh37 coordinates were lifted uniquely to GRCh38 and queried against two
prespecified GTEx v8 brain gene-expression datasets (cortex and DLPFC). All 52 distinct
query responses are cached with hashes. The table reports nominal context only; it does
not assign an eQTL-supported gene or a shared causal signal.

| Region | Pair candidates | Leads | Exact lead eQTL queries / queried | Distinct nominal rows | Lowest nominal p | Gene ID at lowest p |
|---|---:|---|---:|---:|---:|---|
| chr2:57487593-58487593 | 1 | rs11682175 | 2/2 | 13 | 0.149638 | ENSG00000273063 |
| chr2:135040546-136040546 | 1 | rs7599054 | 2/2 | 18 | 0.000313468 | ENSG00000152128 |
| chr3:52219816-53219816 | 1 | rs10933 | 2/2 | 101 | 0.000423087 | ENSG00000163938 |
| chr3:117142005-118142005 | 1 | rs62264767 | 2/2 | 6 | 0.0723162 | ENSG00000185565 |
| chr4:169728549-170728549 | 1 | rs4692709 | 2/2 | 30 | 0.0533106 | ENSG00000129116 |
| chr5:87185500-88185500 | 1 | rs6452785 | 0/2 | 0 | NOT_OBSERVED_AT_EXACT_LEADS | NOT_OBSERVED_AT_EXACT_LEADS |
| chr5:103447968-104481726 | 2 | rs2431108;rs40465 | 4/4 | 4 | 0.0989162 | ENSG00000283462 |
| chr6:104869510-105869510 | 1 | rs4946651 | 2/2 | 23 | 8.53043e-05 | ENSG00000085382 |
| chr7:104429267-105429267 | 1 | rs6466056 | 2/2 | 32 | 0.00660418 | ENSG00000228393 |
| chr7:113559156-114571035 | 2 | rs1476535;rs2894699 | 4/4 | 8 | 0.164376 | ENSG00000135272 |
| chr11:28076261-29076261 | 1 | rs2582894 | 2/2 | 6 | 0.00158934 | ENSG00000169519 |
| chr11:61081764-62081764 | 1 | rs174560 | 2/2 | 156 | 3.06147e-08 | ENSG00000149485 |
| chr11:112877488-113877488 | 1 | rs7111031 | 2/2 | 45 | 0.0166983 | ENSG00000166741 |
| chr13:53365563-54365563 | 1 | rs9536383 | 2/2 | 27 | 0.00150823 | ENSG00000288768 |
| chr14:41646572-42646572 | 1 | rs2121708 | 0/2 | 0 | NOT_OBSERVED_AT_EXACT_LEADS | NOT_OBSERVED_AT_EXACT_LEADS |
| chr14:98137479-99137479 | 1 | rs7160630 | 2/2 | 4 | 0.228114 | ENSG00000285584 |
| chr15:84700520-85700520 | 2 | rs1051168 | 2/2 | 27 | 2.02398e-05 | ENSG00000177082 |
| chr17:1335482-2335482 | 1 | rs4790841 | 0/2 | 0 | NOT_OBSERVED_AT_EXACT_LEADS | NOT_OBSERVED_AT_EXACT_LEADS |
| chr17:43049526-45343136 | 3 | rs17631676;rs199443;rs199528;rs62054437;rs62073157 | 8/10 | 287 | 7.1024e-50 | ENSG00000176681 |
| chr18:52604019-53604019 | 1 | rs4468713 | 2/2 | 12 | 0.286498 | ENSG00000041353 |

The lowest p values are unadjusted minima over varying numbers of genes, lead variants,
and two tissues. They cannot rank regions or justify gene/tissue claims. An absent exact
lead row does not imply that no QTL exists elsewhere in that region. The full v3 evidence
columns remain in the TSV, including replication scope, signed LD, independent LD, and
explicit unrun fine-mapping, coloc, splicing QTL, regulatory, cell-type, and pathway states.
No biologically supported shared region is formally validated or ranked. LAVA local-rg
and final evidence-tier conclusions remain BLOCKED_LAVA under the unchanged canonical
FAILED_QC_NOT_PROMOTED decision.

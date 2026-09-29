# Five-track exploratory molecular credible-set overlap

**EXPLORATORY. Zero LAVA-confirmed regions and zero final shared-locus tiers.**

The four prespecified GTEx v8 eQTL Catalogue datasets cover brain cortex and DLPFC
gene expression and Leafcutter splicing. All 2,686 pair-specific candidate rows remain
in the coverage file. Of those, 2,662 have exact UKB reference alleles and unique
plus-strand GRCh37-to-GRCh38 maps; 24 lack exact reference alleles. Matching requires
the source rsID, lifted coordinate, and unordered allele pair to agree. All four raw
credible-set files are retained and hashed. Molecular PIP is the QTL model's PIP only.

| Region | Candidates in molecular credible sets | Distinct molecular sets | Highest eQTL PIP | Highest sQTL PIP |
|---|---:|---:|---:|---:|
| chr2:57487593-58487593 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr2:135040546-136040546 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr3:52219816-53219816 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr3:117142005-118142005 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr4:169728549-170728549 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr5:87185500-88185500 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr5:103447968-104481726 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr6:104869510-105869510 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr7:104429267-105429267 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr7:113559156-114571035 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr11:28076261-29076261 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr11:61081764-62081764 | 5 | 1 | 0.199166133140605 | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr11:112877488-113877488 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr13:53365563-54365563 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr14:41646572-42646572 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr14:98137479-99137479 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr15:84700520-85700520 | 3 | 1 | 0.00220713919388604 | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr17:1335482-2335482 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |
| chr17:43049526-45343136 | 1686 | 13 | 0.431615618148346 | 0.0332802036973561 |
| chr18:52604019-53604019 | 0 | 0 | NO_EXACT_CANDIDATE_CS_OVERLAP | NO_EXACT_CANDIDATE_CS_OVERLAP |

These are exact PLACO-candidate overlaps with pre-existing **molecular** credible sets.
A PLACO candidate has no established GWAS PIP here, and an overlapping molecular
credible set is not colocalization or proof of a shared causal variant. The maximum
PIP values are descriptive and must not rank shared loci or genes. No absence of QTL
is inferred from an unmatched candidate. Trait fine-mapping, trait-trait coloc,
trait-QTL coloc, independent locus-specific replication, and LAVA local-rg remain
unestablished; final region tiers stay BLOCKED_LAVA.

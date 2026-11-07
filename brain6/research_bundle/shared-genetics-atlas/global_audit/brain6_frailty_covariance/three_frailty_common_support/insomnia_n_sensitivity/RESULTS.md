# Three-frailty comparison: insomnia N sensitivity

Using the raw insomnia per-variant N in place of the archived constant changes the original insomnia–frailty genetic correlations by at most **0.0001 at displayed precision** on the same 1,099,360 SNPs:

| Frailty phenotype | Original r_g | Alternative-N r_g (SE) | Alternative p × 3 |
|---|---:|---:|---:|
| FI | 0.6404 | 0.6405 (0.0236) | 6.10×10⁻¹⁶² |
| FFS | 0.4840 | 0.4839 (0.0270) | 1.45×10⁻⁷¹ |
| HFRS | 0.4355 | 0.4355 (0.0255) | 6.01×10⁻⁶⁵ |

The ranking and three post-hoc paired contrasts are also essentially unchanged. FI minus FFS is 0.1566 (paired SE 0.0281; three-contrast Bonferroni p=7.56×10⁻⁸); HFRS minus FI is −0.2050 (SE 0.0328; adjusted p=1.24×10⁻⁹); HFRS minus FFS is −0.0484 (SE 0.0335; adjusted p=0.445). The two N conditions reuse the same data; this is not independent replication or a formal test of the difference between N conventions. The two correction families each cover only three post-hoc tests, not the broader atlas; nonsignificance does not establish equivalence.

All common-support insomnia rows matched the earlier raw-N map by rsID and allele pair; SNP order, alleles and signed Z were retained. All three new fits retained the same marker count and the insomnia h² delete vectors agreed across pairs. The paired SEs were reconstructed from 200 delete blocks and independently checked in R. These checks do not establish the correct N definition: raw N may be total rather than effective, and no disease-scale h² interpretation follows. FI includes sleep and fatigue; FFS includes exhaustion; HFRS is a different, Finnish hospital-based phenotype with unresolved reference and ascertainment concerns. The N sensitivity does not establish causality, component-excluded genetics, independent cohorts, a local correlation or resolution of the failed LAVA gates.

[Protocol](PROTOCOL.md), [input construction](build.py), [source and matching hashes](build_receipt.json), [LDSC commands](run.sh), [pairwise results](correlations.tsv), [paired contrasts](contrasts.tsv), [calculation](paired.py), [independent check](check.R), and [check log](check.log) document this analysis. Full LDSC logs and delete vectors are in this folder. The inputs and analyses were selected after earlier results, and their agreement is a bounded numerical sensitivity rather than proof that either N choice is correct.

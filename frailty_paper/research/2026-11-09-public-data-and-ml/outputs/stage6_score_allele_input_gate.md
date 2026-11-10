# Score input gate: multiallelic reference sites

A fresh scan of the two local score files on 29 October 2026 finds no A/T or C/G *listed two-allele contrasts* in either prototype. This does not establish the original strand or genome build. Reference-allele listings are those already recorded in the score files, not new genotype-level evidence.

- 15-marker weighted prototype: 9 of 15 sites have more than two alleles in the stored GRCh37 Ensembl reference listing: rs113851554, rs8067974, rs1909650, rs11152363, rs323509, rs4886140, rs35881094, rs11790060 and rs10865954.
- 22-marker count prototype: 14 of 22 sites have more than two alleles in its stored reference listing: rs28528997, rs2815764, rs2155292, rs817509, rs1909650, rs12671330, rs6452792, rs7691487, rs4886140, rs2604551, rs6058750, rs17312861, rs7731706 and rs11636381.

These counts flag reference-site complexity, not observed third-allele frequency in a target cohort. Before scoring, match genome build, chromosome, position, full reference and alternate alleles, the exact source effect contrast and target representation; inspect split multiallelic records and genotype or dosage handling. A same-rsID or same-position match alone cannot establish equivalence. In particular, counting one allele without checking how a third allele was encoded in the discovery model can silently change the comparison group. Prespecify exclusion or a separately justified mapping when the contrast cannot be verified; do not impute it from the allele label alone.

This audit does not resolve the inferred FFS effect direction or demonstrate predictive validity. Source files: `unvalidated_15_variant_score_prototype.csv` and `unvalidated_22_shared_direction_score.csv`; reference annotations and source consistency are described in `candidate_allele_consequence_review.md` and `ffs_effect_allele_source_consistency_audit.md`.

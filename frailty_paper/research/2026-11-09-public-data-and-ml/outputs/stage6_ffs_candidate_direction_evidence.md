# FFS effect direction: candidate-specific coverage check

The prior [37-locus COJO comparison](ffs_effect_allele_source_consistency_audit.md) supports the Ye source file's A1-as-effect-allele convention: all 37 published lead records matched the A1 allele and beta to the reported COJO refA and marginal effect. The [GCTA documentation](https://yanglab.westlake.edu.cn/software/gcta) defines the relevant COJO fields; the [FFS source study](https://www.kp4cd.org/index.php/node/1238) is a UK Biobank genome-wide analysis.

A direct rsID intersection of those 37 audited lead identifiers with the current project score files finds **zero overlaps among the 15 weighted-prototype markers and zero among the 22 shared-count markers**. Thus, the audit is evidence about file-wide coding consistency, not a direct author-table confirmation at any selected score marker. This check is an identifier comparison, not a locus/LD comparison; it does not exclude nearby lead variants or other ways of matching a record.

Maintain the inferred-direction qualifier for individual FFS marker weights. Exact source contrast, build, multiallelic representation and effect orientation must be resolved in the intended target dataset. The 22-marker count does not use FFS beta as a score weight, although FFS direction is relevant to its descriptive interpretation. FFS association coding and the separate TwinGene FI sign problem are distinct issues. The 37 lead checks and PGS score audit do not establish independent score prediction or participant non-overlap.

Reproduce from the three current CSV files using [the intersection script](../work/ffs_candidate_direction_coverage.py). No source or score was changed.

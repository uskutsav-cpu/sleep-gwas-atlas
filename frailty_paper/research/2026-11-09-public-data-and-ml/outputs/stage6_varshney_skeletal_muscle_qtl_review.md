# Skeletal-muscle multiome: candidate molecular follow-up

The [Varshney et al. skeletal-muscle multiome study](https://doi.org/10.1038/s41588-026-02786-6) provides cell-resolved expression and chromatin QTL fine-mapping data from muscle biopsies. Its [public Zenodo release](https://zenodo.org/records/19670721) is CC BY 4.0. This is molecular evidence, not independent testing of the project's frailty or insomnia scores.

We scanned the full released molecular credible-set archives against the 33 unique rsIDs in the 15- and 22-marker lists using previously audited GRCh37→GRCh38 positions and the project allele pairs. In **6,021 expression-QTL credible-set files**, one marker, rs10865954 (GRCh38 `chr3:49174556:T:C`), matched three sets:

| Cell cluster | Gene | Molecular credible-set PIP |
|---|---|---:|
| Endothelial | CCDC36 | 0.0966 |
| Type 1 muscle fiber | USP19 | 0.0610 |
| Type 2x muscle fiber | CCDC36 | 0.00158 |

These are conditional fine-mapping probabilities for molecular signals, not probabilities that the marker causes insomnia or frailty. The source's per-gene permutation tables identify the top eQTL for endothelial CCDC36 and type 1 USP19 as rs11130182, **not** rs10865954. Consequently the source's top-variant beta should not be assigned to the candidate. No expression direction for rs10865954 was established in this screen.

In **75,304 chromatin-QTL credible-set files**, eight candidates matched 58 rows with matching position and allele pair. Their maximum source PIPs across reported molecular signals ranged from 0.0103 to 0.181. These include rs7734045 in eight sets (FAP and fiber clusters; max PIP 0.0115) and rs12536335 in two type 1 fiber sets (max PIP 0.0516). Candidate rs12524473's highest molecular PIP was 0.181, but its extreme BBJ insomnia estimate remains unreliable because the variant was nearly fixed in that study. A PIP is relative to each set's modeled alternatives, and multiple peaks or clusters may share the same underlying signal; rows are not independent replications.

The source-provided GWAS–molecular colocalization tables contain no phenotype labelled insomnia or frailty. Thus these hits are not colocalization with our insomnia or FI/FFS GWAS, evidence for a causal gene, or evidence of cell-specific disease effects. Promoter/peak identity, gene linkage, per-cell sample sizes, molecular and disease study overlap, variant-level QTL effects, matched locus coverage, LD and priors must be checked before a colocalization analysis. The source README points to [CMDGA](https://cmdga.org/) for full cis-QTL statistics, which have not yet been obtained or harmonized for this analysis. The study's molecular credible sets are filtered: no match does not prove a variant has no association, no expression, or no relevant muscle effect.

The expression and chromatin archives were verified against the Zenodo MD5 checksums (`d92556bf8a102acb5fe174e501864a0b` and `eb0a644d5b2f0e46075ae9e8bf4e7d1d`); the smaller downloaded colocalization, expression and permutation files also matched their advertised MD5s. The complete [expression screen receipt](stage6_varshney_muscle_candidate_screen.json) and [chromatin screen receipt](stage6_varshney_muscle_caqtl_candidate_screen.json) include individual set paths, variant identities, PIPs and coverage. Reproduce with `python3 work/varshney_candidate_screen.py` and `python3 work/varshney_caqtl_candidate_screen.py`; the reference mappings and archives must be present. Exact coordinate/pair agreement is conditional on the project build and does not resolve all multiallelic representations.

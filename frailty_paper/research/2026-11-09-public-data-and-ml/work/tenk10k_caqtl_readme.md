This subfolder contains the full caQTL summary statistics for common variants using TensorQTL
Both raw and corrected files are provided

Raw: TenK10K.cis_qtl_pairs.chr*.csv

Corrected: TenK10K.sig_cis_qtl_pairs.chr*.csv

The format for raw association results is as follows:

```
phenotype_id	variant_id	start_distance	af	ma_samples	ma_count	pval_nominal	slopeslope_se
chr21:5065403-5066095	21:5031125:CT:C	-34625	0.14642082	220	270	2.6139886192504902e-11	0.37603652	0.055649046
```

Column interpretation
```
phenotype_id: peak ID, containing chr number and genomic interval
variant_id: variant ID, coded by chr, base pair position, reference allele, and alternative allele
start_distance: distance between variant and peak center
af: allele frequency of the alternative allele
ma_samples: Number of samples that carry at least one copy of the minor allele
ma_count: Total count of the minor allele across all samples
pval_nominal: Nominal (unadjusted) p-value for the association between genotype and phenotype
slope: Estimated effect size of the variant on the phenotype
slope_se: Standard error of the slope estimate
```

The format for the corrected association results is as follows
```
phenotype_id	num_var	beta_shape1	beta_shape2	true_df	pval_true_df	variant_id	start_distance	end_distance	ma_samples	ma_count	af	pval_nominal	slope	slope_se	pval_perm	pval_beta	qval	pval_nominal_threshold
chr2:10183-10900	3554	1.0192658	243.64749	709.8417	0.014203517676235738	2:29350:G:A	18808	18808	117	122	0.06616052	0.007391494668954611	-0.25270417	0.09411177	0.9647035296470353	0.9681353360666333	0.9998751910963442	2.764277935018445e-05
```

Column interpretation
```
phenotype_id: peak ID, containing chr number and genomic interval
num_var: Number of variants tested in cis for this phenotype (within the specified window).
beta_shape1: First shape parameter (α) of the beta distribution used to approximate the permutation p-value distribution.
beta_shape2: Second shape parameter (β) of the beta distribution used for the same approximation.
true_df: Effective degrees of freedom estimated from the data (accounts for covariates and sample size).
pval_true_df: Nominal p-value adjusted using the estimated degrees of freedom instead of the default model assumption.
variant_id: ID of the top associated variant (lead variant) for this peak.
start_distance: Distance (in base pairs) from the peak center to the lead variant.
end_distance: Distance (in base pairs) from the phenotype end to the lead variant.
ma_samples: Number of samples carrying at least one copy of the minor allele.
ma_count: Total count of the minor allele across all samples.
af: Allele frequency of the minor allele (ma_count divided by total alleles).
pval_nominal: Nominal (unadjusted) p-value from the regression for the lead variant.
slope: Effect size (beta); change in phenotype per additional minor allele.
slope_se: Standard error of the effect size.
pval_perm: Empirical p-value from permutations for the best variant (fraction of permutations with equal or stronger association).
pval_beta: Beta-approximated p-value derived from the permutation distribution (more stable than raw permutation p-values).
qval: FDR-adjusted p-value (across phenotypes); used to determine significance (e.g., q < 0.05).
pval_nominal_threshold: Nominal p-value threshold corresponding to the chosen FDR level for this phenotype (used to declare significant variant–phenotype pairs in nominal scans).
```

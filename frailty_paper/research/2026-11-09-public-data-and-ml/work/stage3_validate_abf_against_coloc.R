source('work/coloc_5.2.3_claudia.R')
tests <- list(
 list(path='work/stage3_smim32_FinnGen_insomnia_aligned.csv', result='outputs/stage3_smim32_finngen_coloc_abf_sensitivity.csv', name='SMIM32_FinnGen', type='cc', sdy=1, filter='finn'),
 list(path='work/stage3_smim32_Mak_FI_with_n.csv', result='outputs/stage3_smim32_mak_coloc_abf_sensitivity.csv', name='SMIM32_FI', type='quant', sdy=7.336130485576222, filter='Mak_FI'),
 list(path='work/stage3_smim32_Mak_FP_with_n.csv', result='outputs/stage3_smim32_mak_coloc_abf_sensitivity.csv', name='SMIM32_FP', type='quant', sdy=0.8216409060508033, filter='Mak_FP'),
 list(path='work/stage3_smad5_FinnGen_insomnia_aligned.csv', result='outputs/stage3_smad5_splice_coloc_abf_sensitivity.csv', name='SMAD5_FinnGen', type='cc', sdy=1, filter='FinnGen_insomnia'),
 list(path='work/stage3_smad5_Mak_FI_aligned.csv', result='outputs/stage3_smad5_splice_coloc_abf_sensitivity.csv', name='SMAD5_FI', type='quant', sdy=7.336130485576222, filter='Mak_FI'),
 list(path='work/stage3_smad5_Mak_FP_aligned.csv', result='outputs/stage3_smad5_splice_coloc_abf_sensitivity.csv', name='SMAD5_FP', type='quant', sdy=0.8216409060508033, filter='Mak_FP'),
 list(path='work/stage3_smad5_FinnGen_insomnia_maf005_aligned.csv', result='outputs/stage3_smad5_splice_finngen_maf005_sensitivity.csv', name='SMAD5_FinnGen_MAF005', type='cc', sdy=1, filter='FinnGen_insomnia'))
out <- list()
for (t in tests) {
 d <- read.csv(t$path); q <- approx.bf.estimates(d$qtl_beta_alt/d$qtl_se,d$qtl_se^2,'quant',sdY=1)$lABF
 g <- approx.bf.estimates(d$gwas_beta_alt/d$gwas_se,d$gwas_se^2,t$type,sdY=t$sdy)$lABF
 x <- unname(combine.abf(q,g,1e-4,1e-4,1e-5,quiet=TRUE)); y <- read.csv(t$result)
 if ('trait' %in% colnames(y)) y <- y[y$trait==t$filter,]
 if ('gwas_prior_scale_factor' %in% colnames(y)) y <- y[y$gwas_prior_scale_factor==1,]
 prior_col <- if ('qtl_prior_effect_sd' %in% colnames(y)) 'qtl_prior_effect_sd' else 'qtl_effect_prior_sd'
 y <- y[y$p12==1e-5 & y[[prior_col]]==0.15,]; stopifnot(nrow(y)==1)
 delta <- max(abs(x-as.numeric(y[1,c('H0','H1','H2','H3','H4')])));stopifnot(delta<1e-10)
 out[[length(out)+1]] <- data.frame(analysis=t$name,n=nrow(d),max_absolute_pp_difference=delta)
}
a<-do.call(rbind,out);write.csv(a,'outputs/stage3_coloc_reference_formula_check.csv',row.names=FALSE);print(a)

source('work/coloc_5.2.3_claudia.R')
tests <- list(list(name='Mak_FI',type='quant',sdy=7.336130485576222),list(name='Mak_FP',type='quant',sdy=0.8216409060508033),list(name='FinnGen_insomnia',type='cc',sdy=1))
y<-read.csv('outputs/stage4_cd4_smad5_coloc_abf_sensitivity.csv');out<-list()
for (t in tests) {
 d<-read.csv(paste0('work/stage4_cd4_smad5_',t$name,'_qc.csv'))
 q<-approx.bf.estimates(d$qtl_beta_alt/d$qtl_se,d$qtl_se^2,'quant',sdY=1)$lABF
 g<-approx.bf.estimates(d$gwas_beta_alt/d$gwas_se,d$gwas_se^2,t$type,sdY=t$sdy)$lABF
 a<-unname(combine.abf(q,g,1e-4,1e-4,1e-5,quiet=TRUE)); b<-y[y$trait==t$name & y$p12==1e-5 & y$qtl_prior_sd==.15,];stopifnot(nrow(b)==1)
 delta<-max(abs(a-as.numeric(b[1,c('H0','H1','H2','H3','H4')])));stopifnot(delta<1e-10)
 out[[length(out)+1]]<-data.frame(trait=t$name,n=nrow(d),max_abs_difference=delta)
}
a<-do.call(rbind,out);write.csv(a,'outputs/stage4_cd4_smad5_coloc_formula_check.csv',row.names=FALSE);print(a)

source('work/coloc_5.2.3_claudia.R')
y<-read.csv('outputs/stage3_original_gwas_coloc_sensitivity.csv');out<-list()
for(an in unique(y$analysis)) for(trait in c('Atkins_FI','Jansen_insomnia')){
 d<-read.csv(paste0('work/stage3_',an,'_',trait,'_aligned.csv'));z<-y[y$analysis==an & y$trait==trait & y$p12==1e-5 & y$prior_scale==1,];stopifnot(nrow(z)==1)
 qt<-approx.bf.estimates(d$qtl_beta_alt/d$qtl_se,d$qtl_se^2,'quant',sdY=1)$lABF
 type<-if(trait=='Atkins_FI') 'quant' else 'cc';sdy<-if(type=='quant') z$gwas_prior_sd/.15 else 1
 gw<-approx.bf.estimates(d$gwas_beta_alt/d$gwas_se,d$gwas_se^2,type,sdY=sdy)$lABF
 x<-unname(combine.abf(qt,gw,1e-4,1e-4,1e-5,quiet=TRUE));delta<-max(abs(x-as.numeric(z[1,c('H0','H1','H2','H3','H4')])));stopifnot(delta<1e-10)
 out[[length(out)+1]]<-data.frame(analysis=an,trait=trait,n=nrow(d),max_abs_difference=delta)
}
a<-do.call(rbind,out);write.csv(a,'outputs/stage3_original_gwas_coloc_formula_check.csv',row.names=FALSE);print(a)

source('work/coloc_5.2.3_claudia.R')
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)<=1,length(args)==0 || args[1]=='--tad');suffix<-if(length(args)) '_tad' else ''
tests<-list(list(name='Mak_FI',type='quant'),list(name='Mak_FP',type='quant'),list(name='FinnGen_insomnia',type='cc'))
y<-read.csv(paste0('outputs/stage6_niagads_smad5_exc_abf_sensitivity',suffix,'.csv')); out<-list()
for(t in tests){
 d<-read.csv(paste0('outputs/stage6_niagads_smad5_exc_',t$name,'_aligned',suffix,'.csv'))
 d<-d[abs(d$qtl_alt_af-d$gwas_af_alt)<=.1,]
 b<-y[y$trait==t$name & y$filter=='abs_af_difference_le_0.1' & y$p12==1e-5 & y$qtl_prior_sd==.15,]
 stopifnot(nrow(b)==1,nrow(d)==b$n)
 sdy<-if(t$type=='cc') 1 else b$gwas_prior_sd/.15
 q<-approx.bf.estimates(d$qtl_beta_alt/d$qtl_se,d$qtl_se^2,'quant',sdY=1)$lABF
 g<-approx.bf.estimates(d$gwas_beta_alt/d$gwas_se,d$gwas_se^2,t$type,sdY=sdy)$lABF
 h<-unname(combine.abf(q,g,1e-4,1e-4,1e-5,quiet=TRUE))
 delta<-max(abs(h-as.numeric(b[1,c('H0','H1','H2','H3','H4')]))); stopifnot(delta<1e-10)
 out[[length(out)+1]]<-data.frame(trait=t$name,n=nrow(d),max_abs_difference=delta)
}
a<-do.call(rbind,out);write.csv(a,paste0('outputs/stage6_niagads_smad5_exc_formula_check',suffix,'.csv'),row.names=FALSE);print(a)

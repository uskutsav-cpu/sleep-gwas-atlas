source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p <- 'outputs/shared-genetics-atlas/variant_context/chr5_dense_qtl'
out<-list();tops<-list();i<-0L
for(tissue in c('Testis','Brain_Frontal_Cortex_BA9')) for(gene in c('ENSG00000251574')) for(trait in c('insomnia','adhd')) {
 d<-read.delim(file.path(p,paste0(tissue,'_',gene,'_',trait,'_matched.tsv')),stringsAsFactors=FALSE)
 if(nrow(d)<500){i<-i+1L;out[[i]]<-data.frame(tissue=tissue,gene=gene,trait=trait,n=nrow(d),prior='NOT_ESTIMABLE',p12=NA,H0=NA,H1=NA,H2=NA,H3=NA,H4=NA);next}
 stopifnot(all(is.finite(d$qtl_beta)),all(is.finite(d$qtl_se)),all(d$qtl_se>0),all(is.finite(d$gwas_beta_alt)),all(is.finite(d$gwas_se)),all(d$gwas_se>0),!anyDuplicated(d$variant))
 for(prior in c('default','narrow','wide')){
  if(prior=='default'){
   l1<-approx.bf.estimates(d$qtl_beta/d$qtl_se,d$qtl_se^2,'quant')$lABF
   l2<-approx.bf.estimates(d$gwas_beta_alt/d$gwas_se,d$gwas_se^2,'cc')$lABF
  }else{
   sd<-if(prior=='narrow') .1 else .3;calc<-function(b,se){r<-sd^2/(sd^2+se^2);.5*(log1p(-r)+r*(b/se)^2)}
   l1<-calc(d$qtl_beta,d$qtl_se);l2<-calc(d$gwas_beta_alt,d$gwas_se)
  }
  for(p12 in c(1e-5,1e-6)){
   z<-combine.abf(l1,l2,1e-4,1e-4,p12,quiet=TRUE);stopifnot(abs(sum(z)-1)<1e-12)
   i<-i+1L;out[[i]]<-data.frame(tissue=tissue,gene=gene,trait=trait,n=nrow(d),prior=prior,p12=p12,H0=z[1],H1=z[2],H2=z[3],H3=z[4],H4=z[5],row.names=NULL)
  }
  if(prior=='default'){d$conditional_h4<-exp(l1+l2-logsum(l1+l2));d<-d[order(-d$conditional_h4),];d$credible95<-cumsum(d$conditional_h4)-d$conditional_h4<.95;d$tissue<-tissue;d$gene<-gene;d$trait<-trait;tops[[length(tops)+1L]]<-d}
 }
}
write.table(do.call(rbind,out),file.path(p,'posteriors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
if(length(tops))write.table(do.call(rbind,tops),file.path(p,'conditional_h4.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(do.call(rbind,out))

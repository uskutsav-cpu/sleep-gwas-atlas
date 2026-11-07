# Exploratory, post-selected one-signal ABF; input source P values follow a
# t-distribution, so beta/SE Gaussian ABFs are an approximation.
source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p <- 'outputs/shared-genetics-atlas/variant_context/siraj_2026_mpra/chr5_brain_xqtl/unfiltered_readiness'
data <- '/Volumes/Brain6Work/xqtl-chr5/methylation_candidate'
results<-list();scores<-list()
for(trait in c('insomnia','adhd')){
 d<-read.delim(file.path(data,paste0(trait,'.tsv')),stringsAsFactors=FALSE)
 stopifnot(nrow(d)>500,all(d$qtl_se>0),all(d$gwas_se>0),!anyDuplicated(d$variant))
 for(sd in c(.001,.005,.01,.02)){
  v<-d$qtl_se^2;r<-sd^2/(sd^2+v);l1<-.5*(log1p(-r)+r*(d$qtl_beta/d$qtl_se)^2)
  l2<-approx.bf.estimates(d$gwas_beta_alt/d$gwas_se,d$gwas_se^2,'cc')$lABF
  for(p12 in c(1e-5,1e-6)){
   post<-combine.abf(l1,l2,1e-4,1e-4,p12,quiet=TRUE)
   results[[length(results)+1]]<-data.frame(trait=trait,n=nrow(d),qtl_prior_sd=sd,p12=p12,H0=post[1],H1=post[2],H2=post[3],H3=post[4],H4=post[5],row.names=NULL)
  }
  conditional<-exp(l1+l2-logsum(l1+l2));ix<-order(-conditional)[1:10]
  scores[[length(scores)+1]]<-data.frame(trait=trait,qtl_prior_sd=sd,variant=d$variant[ix],rsid=d$gwas_rsid[ix],conditional_h4=conditional[ix])
 }
}
write.table(do.call(rbind,results),file.path(p,'coloc_posteriors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,scores),file.path(p,'conditional_top.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(p,'session_info.txt'))
print(do.call(rbind,results))

p <- 'outputs/shared-genetics-atlas/component_exclusion/schoeler_physical_function/chr5_coloc'
source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
out<-list();i<-0
for(trait in c('insomnia','adhd')) {
 d<-read.delim(file.path(p,paste0(trait,'_matched.tsv')))
 # Post-hoc fixed effect-size-prior sensitivities; do not infer the outcome SD.
 for(w in c(.05,.10,.15,.30)) {
  r<-w^2/(w^2+d$physical_se^2)
  l1<-.5*(log1p(-r)+r*(d$physical_beta/d$physical_se)^2)
  l2<-approx.bf.estimates(d$trait_beta/d$trait_se,d$trait_se^2,'cc')$lABF
  for(p12 in c(1e-5,1e-6)) {
   z<-combine.abf(l1,l2,1e-4,1e-4,p12,quiet=TRUE)
   i<-i+1;out[[i]]<-data.frame(trait=trait,physical_prior_sd=w,p12=p12,H3=z[4],H4=z[5],row.names=NULL)
  }
 }
}
result<-do.call(rbind,out);write.table(result,file.path(p,'prior_sensitivity.tsv'),sep='\t',row.names=FALSE,quote=FALSE);print(result)

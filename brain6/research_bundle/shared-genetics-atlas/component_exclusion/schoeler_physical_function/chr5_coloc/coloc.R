source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p <- 'outputs/shared-genetics-atlas/component_exclusion/schoeler_physical_function/chr5_coloc'
results<-list(); i<-0
for(trait in c('insomnia','adhd')) {
 d<-read.delim(file.path(p,paste0(trait,'_matched.tsv')))
 stopifnot(nrow(d)>=500, !anyDuplicated(d$pos38), all(d$trait_se>0), all(d$physical_se>0))
 l1<-approx.bf.estimates(d$physical_beta/d$physical_se,d$physical_se^2,'quant')$lABF
 l2<-approx.bf.estimates(d$trait_beta/d$trait_se,d$trait_se^2,'cc')$lABF
 for(p12 in c(1e-5,1e-6)) {
  z<-combine.abf(l1,l2,1e-4,1e-4,p12,quiet=TRUE);stopifnot(abs(sum(z)-1)<1e-12)
  i<-i+1;results[[i]]<-data.frame(trait=trait,n=nrow(d),p12=p12,H0=z[1],H1=z[2],H2=z[3],H3=z[4],H4=z[5],min_physical_p=min(d$physical_p),min_trait_p=min(pchisq((d$trait_beta/d$trait_se)^2,1,lower.tail=FALSE)),row.names=NULL)
 }
 d$conditional_h4<-exp(l1+l2-logsum(l1+l2));d<-d[order(-d$conditional_h4),];d$credible95<-cumsum(d$conditional_h4)-d$conditional_h4<.95
 write.table(d,file.path(p,paste0(trait,'_conditional.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
}
result<-do.call(rbind,results);write.table(result,file.path(p,'posteriors.tsv'),sep='\t',row.names=FALSE,quote=FALSE);print(result)

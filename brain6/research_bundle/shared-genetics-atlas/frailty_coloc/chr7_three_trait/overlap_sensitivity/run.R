.libPaths(c('/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
library(hyprcoloc)
p<-'outputs/shared-genetics-atlas/frailty_coloc/chr7_three_trait/overlap_sensitivity'
x<-read.delim('/Volumes/Brain6Work/chr7-three-trait-input.tsv')
stopifnot(nrow(x)==1055, !anyDuplicated(x$rsid))
b<-cbind(FI=x$fi_beta,insomnia=x$insomnia_beta,ADHD=x$adhd_beta)
s<-cbind(FI=x$fi_se,insomnia=x$insomnia_se,ADHD=x$adhd_se)
results<-list()
for (rho in c(0,-.3,-.1,-.001,.001,.1,.3,.5)) {
  mat<-diag(3); mat[1,2]<-mat[2,1]<-rho
  stopifnot(min(eigen(mat,symmetric=TRUE,only.values=TRUE)$values)>0)
  cat('rho',rho,'\n'); flush.console()
  result<-tryCatch(hyprcoloc(b,s,binary.outcomes=c(0,1,1),trait.names=colnames(b),snp.id=x$rsid,prior.1=1e-4,prior.c=.02,reg.thresh=.5,align.thresh=.5,ind.traits=(rho==0),trait.cor=mat,sample.overlap=matrix(1,3,3),bb.alg=TRUE,snpscores=TRUE),error=function(e) list(error=conditionMessage(e)))
  results[[as.character(rho)]]<-result
  print(result$results); if(!is.null(result$error)) cat('ERROR',result$error,'\n'); flush.console()
  saveRDS(results,file.path(p,'sensitivity_results.rds'))
}
writeLines(capture.output(sessionInfo()),file.path(p,'session_info.txt'))

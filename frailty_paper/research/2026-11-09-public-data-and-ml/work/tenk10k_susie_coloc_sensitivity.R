.libPaths(c('work/stage5_Rlib',.libPaths()));library(coloc)
a<-readRDS('work/stage6_tenk10k_susie_diagnostic/fits.rds');d<-read.csv('work/stage6_tenk10k_susie_diagnostic/aligned.csv');stopifnot(a$qtl$converged,a$fi$converged)
for(k in names(a)) for(m in c('alpha','lbf_variable')) colnames(a[[k]][[m]])<-d$variant_GRCh38
res<-list()
for(p12 in c(1e-6,1e-5,5e-5)){
 x<-coloc.susie(a$qtl,a$fi,p1=1e-4,p2=1e-4,p12=p12)$summary
 stopifnot(nrow(x)==1);x$p12<-p12;res[[length(res)+1]]<-as.data.frame(x)
}
y<-do.call(rbind,res);write.csv(y,'outputs/stage6_tenk10k_susie_coloc_sensitivity.csv',row.names=FALSE);print(y)

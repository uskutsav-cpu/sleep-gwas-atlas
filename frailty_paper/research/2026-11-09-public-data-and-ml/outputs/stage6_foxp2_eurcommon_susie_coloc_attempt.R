.libPaths(c('work/stage5_Rlib',.libPaths()));library(coloc);library(susieR)
h<-readRDS('work/stage6_foxp2_eurcommon_FinnGen_HFRS_susie.rds');res<-list();n<-0
for (nchoice in c('effN','totalN')){
 m<-readRDS(paste0('work/stage6_foxp2_eurcommon_MVP_insomnia_',nchoice,'_susie.rds'))
 for(p12 in c(1e-6,1e-5,5e-5)){
  r<-coloc.susie(m,h,p1=1e-4,p2=1e-4,p12=p12); s<-as.data.frame(r$summary); s$mvp_n_assumption<-nchoice;s$p1<-1e-4;s$p2<-1e-4;s$p12<-p12;n<-n+1;res[[n]]<-s
 }
}
out<-do.call(rbind,res);write.csv(out,'outputs/stage6_foxp2_eurcommon_susie_coloc_prior_sensitivity.csv',row.names=FALSE);print(out)

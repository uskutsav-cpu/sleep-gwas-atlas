.libPaths(c('work/stage5_Rlib',.libPaths()));library(susieR)
files<-c('eurcommon_MVP_insomnia_effN','eurcommon_FinnGen_HFRS','fin_MVP_insomnia_effN','fin_FinnGen_HFRS');out<-list()
for(name in files){f<-readRDS(paste0('work/stage6_foxp2_',name,'_susie.rds'));if(length(f$sets$cs)>0)for(i in seq_along(f$sets$cs)){ids<-f$sets$cs[[i]];rs<-colnames(f$alpha)[ids];out[[length(out)+1]]<-data.frame(model=name,component=f$sets$cs_index[i],rsid=rs,pip=f$pip[ids])}}
write.csv(do.call(rbind,out),'outputs/stage6_foxp2_credible_set_members.csv',row.names=FALSE)

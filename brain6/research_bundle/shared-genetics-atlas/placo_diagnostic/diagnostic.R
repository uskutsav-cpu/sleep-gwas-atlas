args <- commandArgs(trailingOnly=TRUE)
source('/Volumes/Extreme SSD/brain6-work/PLACO-v0.2.0/PLACO_v0.2.0.R')
d <- read.delim(args[1],check.names=FALSE)
# These values are independently calculated from all aligned HM3 support in diagnostic_parameters.json.
settings <- list(FI=c(1.271851342700776,1.4045250827610078,0.08081972497537657), FFS=c(1.2714481467007943,1.5006275551853958,0.052417344112559244))
start<-proc.time()[3]; out<-list()
for (trait in names(settings)) {
 z2<-d[[paste0('Z_',trait)]]; z1<-d$Z_insomnia; par<-settings[[trait]]
 ans<-vapply(seq_len(nrow(d)),function(i) {
  if(z1[i]==0 || z2[i]==0) return(1)
  tryCatch({x<-placo.plus(c(z1[i],z2[i]),VarZ=par[1:2],CorZ=par[3],AbsTol=1e-13)$p.placo.plus; if(length(x)==1 && is.finite(x) && x>=0 && x<=1) x else NA_real_},error=function(e) NA_real_)
 },numeric(1))
 out[[trait]]<-ans; cat(trait,'n=',length(ans),'fail=',sum(is.na(ans)),'p<0.05=',sum(ans<.05,na.rm=TRUE),'elapsed=',proc.time()[3]-start,'\n')
}
d$P_PLACO_FI<-out$FI;d$P_PLACO_FFS<-out$FFS
write.table(d,args[2],quote=FALSE,sep='\t',row.names=FALSE)

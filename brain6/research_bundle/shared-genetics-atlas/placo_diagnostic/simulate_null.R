args<-commandArgs(trailingOnly=TRUE)
suppressMessages(source('/Volumes/Extreme SSD/brain6-work/PLACO-v0.2.0/PLACO_v0.2.0.R'))
settings<-list(FI=c(1.271851342700776,1.4045250827610078,0.08081972497537657),FFS=c(1.2714481467007943,1.5006275551853958,0.052417344112559244))
set.seed(20261027); n<-1000L;result<-list();idx<-0
for (trait in names(settings)) for(scenario in c('both_null','sleep_only','frailty_only')) {
 par<-settings[[trait]];u<-rnorm(n);v<-rnorm(n); z1<-u;z2<-par[3]*u+sqrt(1-par[3]^2)*v
 if(scenario=='sleep_only')z1<-z1+2
 if(scenario=='frailty_only')z2<-z2+2
 p<-vapply(seq_len(n),function(i)tryCatch({x<-placo.plus(c(z1[i],z2[i]),VarZ=par[1:2],CorZ=par[3],AbsTol=1e-13)$p.placo.plus;if(length(x)==1&&is.finite(x)&&x>0&&x<=1)x else NA_real_},error=function(e)NA_real_),numeric(1))
 idx<-idx+1;result[[idx]]<-data.frame(pair=trait,scenario=scenario,Z1=z1,Z2=z2,P=p);cat(trait,scenario,'errors',sum(is.na(p)),'p<.05',sum(p<.05,na.rm=TRUE),'\n');flush.console()
}
write.table(do.call(rbind,result),args[1],sep='\t',quote=FALSE,row.names=FALSE,na='NA')

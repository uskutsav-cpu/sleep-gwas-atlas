args<-commandArgs(trailingOnly=TRUE);suppressMessages(source('/Volumes/Extreme SSD/brain6-work/PLACO-v0.2.0/PLACO_v0.2.0.R'));settings<-list(FI=c(1.271851342700776,1.4045250827610078,0.08081972497537657),FFS=c(1.2714481467007943,1.5006275551853958,0.052417344112559244));set.seed(20261028);n<-1000L;res<-list();k<-0
for(trait in names(settings))for(scenario in c('both_unit','sleep_variance_inflated','frailty_variance_inflated')){
 q<-settings[[trait]];u<-rnorm(n);v<-rnorm(n);z1<-u;z2<-q[3]*u+sqrt(1-q[3]^2)*v;if(scenario=='sleep_variance_inflated')z1<-z1*sqrt(q[1]);if(scenario=='frailty_variance_inflated')z2<-z2*sqrt(q[2]);
 p<-vapply(seq_len(n),function(i)tryCatch({x<-placo.plus(c(z1[i],z2[i]),VarZ=q[1:2],CorZ=q[3],AbsTol=1e-13)$p.placo.plus;if(length(x)==1&&is.finite(x)&&x>0&&x<=1)x else NA_real_},error=function(e)NA_real_),numeric(1));k<-k+1;res[[k]]<-data.frame(pair=trait,scenario=scenario,Z1=z1,Z2=z2,P=p);cat(trait,scenario,'errors',sum(is.na(p)),'p<.05',sum(p<.05,na.rm=TRUE),'\n');flush.console()
}
write.table(do.call(rbind,res),args[1],sep='\t',quote=FALSE,row.names=FALSE,na='NA')

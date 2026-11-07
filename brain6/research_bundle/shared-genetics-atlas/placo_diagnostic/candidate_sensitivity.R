args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1);out<-args[1];d<-read.delim(file.path(out,'candidate_source_matches.tsv'),na.strings=c('','NA'))
source('/Volumes/Extreme SSD/brain6-work/PLACO-v0.2.0/PLACO_v0.2.0.R');par<-list(FI=c(1.271851342700776,1.4045250827610078,0.08081972497537657),FFS=c(1.2714481467007943,1.5006275551853958,0.052417344112559244))
evalp<-function(z,params){if(any(!is.finite(z))||any(z^2>80))return(NA_real_);if(any(z==0))return(1);tryCatch({p<-placo.plus(z,VarZ=params[1:2],CorZ=params[3],AbsTol=1e-13)$p.placo.plus;if(length(p)==1 && is.finite(p) && p>0 && p<=1)p else NA_real_},error=function(e)NA_real_)}
d$p_source_p<-NA_real_;d$p_source_wald<-NA_real_
if(nrow(d))for(i in seq_len(nrow(d))) if(d$sensitivity_eligible[i]=='YES') {d$p_source_p[i]<-evalp(c(d$sleep_z_p[i],d$frailty_z_p[i]),par[[d$pair[i]]]);d$p_source_wald[i]<-evalp(c(d$sleep_z_wald[i],d$frailty_z_wald[i]),par[[d$pair[i]]])}
write.table(d,file.path(out,'candidate_statistic_sensitivity.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA');cat('Total',nrow(d),'eligible',sum(d$sensitivity_eligible=='YES'),'\n')

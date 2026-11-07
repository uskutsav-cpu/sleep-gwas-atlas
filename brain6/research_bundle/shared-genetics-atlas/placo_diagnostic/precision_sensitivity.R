# Directly compare the same official PLACO+ under three source-Z constructions.
args<-commandArgs(trailingOnly=TRUE);limit<-if(length(args)) as.integer(args[1]) else Inf
root<-'outputs/shared-genetics-atlas/placo_diagnostic'; source('/Volumes/Extreme SSD/brain6-work/PLACO-v0.2.0/PLACO_v0.2.0.R')
d<-read.delim(file.path(root,'precision_audit_sample.tsv'));a<-d[d$source=='S',];a<-a[order(a$SNP),];if(is.finite(limit))a<-head(a,limit)
settings<-list(FI=c(1.271851342700776,1.4045250827610078,0.08081972497537657),FFS=c(1.2714481467007943,1.5006275551853958,0.052417344112559244))
out<-list();j<-0L
for(pair in c('FI','FFS')) {
 b<-d[d$source==pair,];x<-merge(a,b,by=c('SNP','A1','A2'),suffixes=c('_s','_f'));stopifnot(nrow(x)==nrow(a),!anyDuplicated(x$SNP));par<-settings[[pair]]
 for(i in seq_len(nrow(x))) {
  pair_eval<-function(z){if(any(!is.finite(z))||any(z^2>80))return(NA_real_);if(any(z==0))return(1);tryCatch({p<-placo.plus(z,VarZ=par[1:2],CorZ=par[3],AbsTol=1e-13)$p.placo.plus;if(length(p)==1 && is.finite(p) && p>0 && p<=1)p else NA_real_},error=function(e) NA_real_)}
  zrounded<-c(x$z_munged_s[i],x$z_munged_f[i]);zp<-c(x$z_from_reported_p_and_beta_sign_s[i],x$z_from_reported_p_and_beta_sign_f[i]);zwald<-c(x$z_wald_beta_se_s[i],x$z_wald_beta_se_f[i]);j<-j+1L
  out[[j]]<-data.frame(pair=pair,SNP=x$SNP[i],z_sleep_rounded=zrounded[1],z_frailty_rounded=zrounded[2],z_sleep_p=zp[1],z_frailty_p=zp[2],z_sleep_wald=zwald[1],z_frailty_wald=zwald[2],p_rounded=pair_eval(zrounded),p_source_p=pair_eval(zp),p_wald=pair_eval(zwald))
 }
 cat(pair,'completed',nrow(x),'\n');flush.console()
}
o<-do.call(rbind,out);fname<-if(is.finite(limit)) paste0('precision_sensitivity_smoke_',limit,'.tsv') else 'precision_sensitivity.tsv';write.table(o,file.path(root,fname),row.names=FALSE,quote=FALSE,sep='\t');for(pair in c('FI','FFS')) {x<-o[o$pair==pair,];cat(pair,'n',nrow(x),'failed',sum(!complete.cases(x)),'median_abs_log10_round_vs_p',median(abs(log10(x$p_rounded)-log10(x$p_source_p)),na.rm=TRUE),'median_abs_log10_p_vs_wald',median(abs(log10(x$p_source_p)-log10(x$p_wald)),na.rm=TRUE),'discordant_0.05',sum((x$p_source_p<.05)!=(x$p_wald<.05),na.rm=TRUE),'\n')}

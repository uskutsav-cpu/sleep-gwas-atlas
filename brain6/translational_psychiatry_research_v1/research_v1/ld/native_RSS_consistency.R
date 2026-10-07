lib<-normalizePath('work/R-scientific-lib');.libPaths(c(lib,.libPaths()))
stopifnot(as.character(packageVersion('susieR'))=='0.14.2',as.character(packageVersion('coloc'))=='5.2.3')
library(susieR)
here<-'brain6/translational_psychiatry_research_v1/research_v1/ld'
local<-'work/ld_genotypes_research_v1/C_finemap_inputs'
source(file.path(here,'native_runtime_adapter.R'))
checks<-read.delim(file.path(here,'finemap_input_validation.tsv'),check.names=FALSE)
if(!all(checks$input_gate_status=='PASS'))stop('Input QC hold')
a<-read.delim(file.path(local,'insomnia_C_shared.tsv'));b<-read.delim(file.path(local,'adhd_C_shared.tsv'))
stopifnot(identical(a$variant_key,b$variant_key),nrow(a)>=500)
p<-nrow(a);R<-matrix(readBin(file.path(local,'shared_signed_LD.float64.bin'),what=numeric(),n=p*p,size=8,endian='little'),nrow=p,ncol=p)
stopifnot(all(is.finite(R)),max(abs(R-t(R)))<=1e-12,max(abs(diag(R)-1))<=1e-12,max(abs(R))<=1+1e-12)
e<-eigen(R,symmetric=TRUE);stopifnot(min(e$values)>=-1e-8*max(e$values));attr(R,'eigen')<-e
outputs<-list();flags<-list()
for(trait in c('insomnia','adhd')) {
 d<-if(trait=='insomnia')a else b;n<-median(d$N_eff);z<-d$Z_ALT
 s<-estimate_s_rss(z,R,n=n,method='null-mle');s_inf<-rss_large_n_consistency(z,R)
 k<-kriging_rss(z,R,n=n,s=s)
 k$conditional_dist<-extract_kriging_table(k,p)
 write.table(k$conditional_dist,file.path(local,paste0(trait,'_kriging_diagnostics_LOCAL_ONLY.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 cat('DIAGNOSTIC_COLUMNS ',paste(names(k$conditional_dist),collapse=','),'\n')
 zz<-k$conditional_dist[['z']];if(is.null(zz))zz<-k$conditional_dist[['zscore']]
 lr<-k$conditional_dist[['logLR']];stopifnot(length(lr)==p,length(zz)==p)
 flag<-is.finite(lr)&lr>2&abs(zz)>2
 result<-data.frame(trait=trait,region='chr5:103447968-104447968',n_shared_snps=p,n_effective_median=n,s_null_mle=s,s_null_mle_N_Inf=s_inf,n_allele_switch_flags=sum(flag),max_logLR=max(lr),consistency_status=if(s<=.10&&!any(flag))'PASS_OPERATIONAL_CONSISTENCY_GATE' else 'HOLD_GWAS_LD_CONSISTENCY',criteria='s<=.10_AND_ZERO_logLRgt2_absZgt2',method='susieR0.14.2_native_estimate_s_rss_and_kriging_rss',independent_replication=FALSE)
 outputs[[trait]]<-result
 if(any(flag))flags[[trait]]<-data.frame(trait=trait,variant_key=d$variant_key[flag],SNP=d$SNP[flag],logLR=lr[flag],absZ_gt2=TRUE,interpretation='POTENTIAL_SOURCE_LD_OR_ALLELE_INCOMPATIBILITY_NOT_CONFIRMED_ALLELE_ERROR')
 saveRDS(list(trait=trait,s=s,s_inf=s_inf,kriging_table=k$conditional_dist,LD_eigen_summary=c(min=min(e$values),max=max(e$values)),input_keys=d$variant_key),file.path(local,paste0(trait,'_RSS_consistency.rds')))
 print(result)
 print(warnings())
}
write.table(do.call(rbind,outputs),file.path(here,'native_RSS_consistency_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
if(length(flags)) {write.table(do.call(rbind,flags),file.path(here,'native_RSS_allele_switch_flags.tsv'),sep='\t',row.names=FALSE,quote=FALSE)} else {writeLines('trait\tvariant_key\tSNP\tlogLR\tabsZ_gt2\tinterpretation',file.path(here,'native_RSS_allele_switch_flags.tsv'))}
writeLines(capture.output(sessionInfo()),file.path(here,'native_RSS_sessionInfo.txt'))

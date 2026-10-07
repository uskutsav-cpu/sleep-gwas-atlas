lib<-normalizePath('work/R-scientific-lib');.libPaths(c(lib,.libPaths()))
library(susieR);library(coloc)
stopifnot(as.character(packageVersion('susieR'))=='0.14.2',as.character(packageVersion('coloc'))=='5.2.3')
here<-'brain6/translational_psychiatry_research_v1/research_v1/ld';local<-'work/ld_genotypes_research_v1/C_finemap_inputs'
checks<-read.delim(file.path(here,'finemap_input_validation.tsv'),check.names=FALSE);consistency<-read.delim(file.path(here,'native_RSS_consistency_summary.tsv'))
stopifnot(all(checks$input_gate_status=='PASS'),all(consistency$consistency_status=='PASS_OPERATIONAL_CONSISTENCY_GATE'),all(consistency$s_null_mle<=.10),all(consistency$n_allele_switch_flags==0))
data<-list(insomnia=read.delim(file.path(local,'insomnia_C_shared.tsv')),adhd=read.delim(file.path(local,'adhd_C_shared.tsv')))
stopifnot(identical(data$insomnia$variant_key,data$adhd$variant_key));keys<-data$insomnia$variant_key;p<-length(keys)
R<-matrix(readBin(file.path(local,'shared_signed_LD.float64.bin'),what=numeric(),n=p*p,size=8,endian='little'),nrow=p,ncol=p,dimnames=list(keys,keys))
stopifnot(p>=500,all(is.finite(R)),max(abs(R-t(R)))<=1e-12,max(abs(diag(R)-1))<=1e-12)
models<-data.frame(model=c('PRIMARY_L10_MEDIAN_NEFF','SENSITIVITY_L1_MEDIAN_NEFF','SENSITIVITY_L5_MEDIAN_NEFF','SENSITIVITY_L10_N_INFINITY'),L=c(10,1,5,10),use_n=c(TRUE,TRUE,TRUE,FALSE))
fit_summaries<-list();pip_rows<-list();cs_rows<-list();elbo_rows<-list();coloc_rows<-list();fit_no<-0;coloc_no<-0
save_tables<-function() {
 if(length(fit_summaries))write.table(do.call(rbind,fit_summaries),file.path(here,'native_finemap_fit_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 if(length(pip_rows))write.table(do.call(rbind,pip_rows),file.path(here,'native_finemap_variant_PIPs.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 if(length(cs_rows))write.table(do.call(rbind,cs_rows),file.path(here,'native_finemap_credible_sets.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 if(length(elbo_rows))write.table(do.call(rbind,elbo_rows),file.path(here,'native_finemap_ELBO_history.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 if(length(coloc_rows))write.table(do.call(rbind,coloc_rows),file.path(here,'native_coloc_all_signal_pairs.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
}
for(mi in seq_len(nrow(models))) {
 model<-models$model[mi];L<-models$L[mi];fits<-list()
 for(trait in names(data)) {
  d<-data[[trait]];z<-setNames(d$Z_ALT,keys);n<-median(d$N_eff);warns<-character()
  fitargs<-list(z=z,R=R,L=L,estimate_residual_variance=FALSE,coverage=.95,min_abs_corr=.5,max_iter=1000,tol=1e-3,refine=FALSE)
  if(models$use_n[mi])fitargs$n<-n
  fit<-withCallingHandlers(do.call(susieR::susie_rss,fitargs),warning=function(w){warns<<-c(warns,conditionMessage(w))})
  colnames(fit$lbf_variable)<-keys;colnames(fit$alpha)<-keys;names(fit$pip)<-keys
  saveRDS(fit,file.path(local,paste0(trait,'_',model,'_native_fit.rds')))
  csn<-if(is.null(fit$sets$cs))0 else length(fit$sets$cs);drop<-if(length(fit$elbo)>1)min(diff(fit$elbo)) else 0
  good<-isTRUE(fit$converged)&&all(is.finite(fit$pip))&&all(fit$pip>=0&fit$pip<=1)&&drop>=-1e-6
  status<-if(!good)'FAIL_CONVERGENCE_OR_NUMERIC' else if(csn==0)'NO_QUALIFYING_CREDIBLE_SET' else 'EXPLORATORY_REFERENCE_QUALIFIED_FIT'
  fit_no<-fit_no+1
  fit_summaries[[fit_no]]<-data.frame(model=model,trait=trait,region='chr5:103447968-104447968',n_snps=p,L=L,n_mode=if(models$use_n[mi])'MEDIAN_NEFF_BINARY_APPROXIMATION' else 'N_INFINITY_Z_ONLY_APPROXIMATION',n_effective_if_supplied=if(models$use_n[mi])n else NA,converged=isTRUE(fit$converged),n_iterations=length(fit$elbo),min_ELBO_change=drop,n_credible_sets=csn,max_PIP=max(fit$pip),status=status,warnings=paste(warns,collapse=' | '),independent_replication=FALSE)
  pip_rows[[fit_no]]<-data.frame(model=model,trait=trait,variant_key=keys,SNP=d$SNP,PIP=fit$pip)
  elbo_rows[[fit_no]]<-data.frame(model=model,trait=trait,iteration=seq_along(fit$elbo),ELBO=fit$elbo)
  if(csn)for(j in seq_len(csn)) {
   ix<-fit$sets$cs[[j]];purity<-fit$sets$purity[j,,drop=FALSE]
   cs_rows[[length(cs_rows)+1]]<-data.frame(model=model,trait=trait,credible_set=names(fit$sets$cs)[j],effect_index=fit$sets$cs_index[j],CS_size=length(ix),purity_min_abs_corr=purity$min.abs.corr,purity_mean_abs_corr=purity$mean.abs.corr,purity_median_abs_corr=purity$median.abs.corr,variant_key=keys[ix],SNP=d$SNP[ix],variant_PIP=fit$pip[ix],within_signal_posterior=fit$alpha[fit$sets$cs_index[j],ix])
  }
  fits[[trait]]<-list(fit=fit,good=good,csn=csn,status=status)
  print(fit_summaries[[fit_no]]);save_tables()
 }
 for(p12 in c(1e-5,1e-6,5e-5)) {
  coloc_no<-coloc_no+1
  if(!all(vapply(fits,function(f)isTRUE(f$good)&&f$csn>0,logical(1)))) {
   coloc_rows[[coloc_no]]<-data.frame(model=model,p1=1e-4,p2=1e-4,p12=p12,status='NOT_ESTIMABLE_FIT_OR_CREDIBLE_SET_HOLD',nsnps=p,hit1=NA,hit2=NA,PP.H0.abf=NA,PP.H1.abf=NA,PP.H2.abf=NA,PP.H3.abf=NA,PP.H4.abf=NA,idx1=NA,idx2=NA,H4_over_H3_plus_H4=NA,descriptive_support=FALSE,independent_replication=FALSE)
  } else {
   co<-coloc::coloc.susie(fits$insomnia$fit,fits$adhd$fit,p1=1e-4,p2=1e-4,p12=p12)
   saveRDS(co,file.path(local,paste0(model,'_p12_',format(p12,scientific=TRUE),'_coloc_susie.rds')))
   ss<-as.data.frame(co$summary);ss$model<-model;ss$p1<-1e-4;ss$p2<-1e-4;ss$p12<-p12;ss$status<-'EXPLORATORY_NATIVE_COLOC_SUSIE';ss$H4_over_H3_plus_H4<-ss$PP.H4.abf/(ss$PP.H3.abf+ss$PP.H4.abf);ss$descriptive_support<-ss$PP.H4.abf>=.8&ss$H4_over_H3_plus_H4>=.8;ss$independent_replication<-FALSE
   coloc_rows[[coloc_no]]<-ss[,c('model','p1','p2','p12','status','nsnps','hit1','hit2','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf','idx1','idx2','H4_over_H3_plus_H4','descriptive_support','independent_replication')]
   print(coloc_rows[[coloc_no]])
  }
  save_tables()
 }
}
# Independent native single-signal implementation on identical fixed SNP universe.
abf_rows<-list()
for(p12 in c(1e-5,1e-6,5e-5)) {
 d1<-list(beta=data$insomnia$beta_ALT,varbeta=data$insomnia$SE^2,snp=keys,type='cc')
 d2<-list(beta=data$adhd$beta_ALT,varbeta=data$adhd$SE^2,snp=keys,type='cc')
 co<-coloc::coloc.abf(d1,d2,p1=1e-4,p2=1e-4,p12=p12)
 saveRDS(co,file.path(local,paste0('same_universe_p12_',format(p12,scientific=TRUE),'_coloc_abf.rds')))
 ss<-as.data.frame(as.list(co$summary));ss$p1<-1e-4;ss$p2<-1e-4;ss$p12<-p12;ss$model<-'SAME_UNIVERSE_SINGLE_SIGNAL_ABF';ss$H4_over_H3_plus_H4<-ss$PP.H4.abf/(ss$PP.H3.abf+ss$PP.H4.abf);ss$independent_replication<-FALSE;abf_rows[[length(abf_rows)+1]]<-ss
}
write.table(do.call(rbind,abf_rows),file.path(here,'same_universe_native_ABF_sensitivity.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(here,'native_finemap_sessionInfo.txt'))
cat('COMPLETED_FIXED_NATIVE_FITS ',fit_no,' COLOC_INVOCATIONS ',coloc_no,' ABF_INVOCATIONS ',length(abf_rows),'\n')

.libPaths(c('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
stopifnot(requireNamespace('susieR',quietly=TRUE),requireNamespace('coloc',quietly=TRUE))
p <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_dense_qtl/multisignal_pilot'
v<-read.delim(file.path(p,'variants.tsv'),stringsAsFactors=FALSE);ns<-read.delim(file.path(p,'sample_sizes.tsv'),stringsAsFactors=FALSE);stopifnot(identical(v$variant,ns$variant),nrow(v)==1548)
R<-as.matrix(read.delim(gzfile(file.path(p,'ld_matrix.tsv.gz')),header=FALSE,check.names=FALSE));dimnames(R)<-list(v$variant,v$variant)
stopifnot(nrow(R)==nrow(v),ncol(R)==nrow(v),all(is.finite(R)),max(abs(R-t(R)))<1e-9,max(abs(diag(R)-1))<1e-9)
source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
base<-list();i<-0L
for(trait in c('insomnia','adhd')){
 l1<-approx.bf.estimates(v$qtl_beta/v$qtl_se,v$qtl_se^2,'quant')$lABF;l2<-approx.bf.estimates(v[[paste0(trait,'_beta')]]/v[[paste0(trait,'_se')]],v[[paste0(trait,'_se')]]^2,'cc')$lABF
 for(p12 in c(1e-5,1e-6)){z<-combine.abf(l1,l2,1e-4,1e-4,p12,quiet=TRUE);i<-i+1L;base[[i]]<-data.frame(trait=trait,n=nrow(v),p12=p12,H0=z[1],H1=z[2],H2=z[3],H3=z[4],H4=z[5],row.names=NULL)}
}
write.table(do.call(rbind,base),file.path(p,'single_signal_same_support.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
settings<-list(qtl=list(z=v$qtl_beta/v$qtl_se,n=round(median(ns$qtl_N_from_an))),insomnia=list(z=v$insomnia_beta/v$insomnia_se,n=round(median(ns$insomnia_N_AUDIT))),adhd=list(z=v$adhd_beta/v$adhd_se,n=round(median(ns$adhd_N_AUDIT))),adhd_minN=list(z=v$adhd_beta/v$adhd_se,n=round(min(ns$adhd_N_AUDIT))))
set.seed(20261027)
fits<-list();reports<-list();credible<-list();diagnostics<-list();i<-0L
for(name in names(settings)){
 d<-settings[[name]];names(d$z)<-v$variant;notes<-character();start<-proc.time()[['elapsed']];cat('starting',name,'N',d$n,'\n');flush.console()
 result<-tryCatch(withCallingHandlers(susieR::susie_rss(z=d$z,R=R,n=d$n,L=10,estimate_residual_variance=FALSE,max_iter=100,tol=1e-3),warning=function(w){notes<<-c(notes,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(result,'error')){i<-i+1L;reports[[i]]<-data.frame(fit=name,N=d$n,status='ERROR',converged=FALSE,iterations=NA,credible_sets=NA,seconds=proc.time()[['elapsed']]-start,warnings=paste(c(notes,conditionMessage(result)),collapse=' | '));next}
 fits[[name]]<-result;saveRDS(result,file.path(p,paste0('fit_',name,'.rds')))
 ncs<-length(result$sets$cs);i<-i+1L;reports[[i]]<-data.frame(fit=name,N=d$n,status='OK',converged=isTRUE(result$converged),iterations=result$niter,credible_sets=ncs,seconds=proc.time()[['elapsed']]-start,warnings=paste(notes,collapse=' | '));cat('finished',name,'converged',isTRUE(result$converged),'sets',ncs,'\n');flush.console()
 if(ncs>0){for(k in seq_len(ncs)){ix<-result$sets$cs[[k]];purity<-result$sets$purity[k,,drop=FALSE];credible[[length(credible)+1L]]<-data.frame(fit=name,credible_set=names(result$sets$cs)[k],size=length(ix),min_abs_corr=purity[['min.abs.corr']],mean_abs_corr=purity[['mean.abs.corr']],variants=paste(v$variant[ix],collapse=','))}}
}
write.table(do.call(rbind,reports),file.path(p,'fit_summary.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
if(length(credible))write.table(do.call(rbind,credible),file.path(p,'credible_sets.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
for(name in c('qtl','insomnia','adhd')){
 d<-settings[[name]];cat('LD diagnostic',name,'\n');flush.console()
 value<-tryCatch(susieR::estimate_s_rss(z=d$z,R=R,n=d$n),error=function(e)NA_real_);diagnostics[[name]]<-data.frame(fit=name,n=d$n,estimated_s=value)
}
write.table(do.call(rbind,diagnostics),file.path(p,'ld_diagnostics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
comp<-list()
for(trait in c('insomnia','adhd','adhd_minN')){
 if(is.null(fits[['qtl']]) || is.null(fits[[trait]]) || !isTRUE(fits[['qtl']]$converged) || !isTRUE(fits[[trait]]$converged)){comp[[length(comp)+1L]]<-data.frame(trait=trait,p12=NA,status='NOT_ESTIMABLE_FIT',detail='failed or not converged');next}
 for(prior in c(1e-5,1e-6)){
  ans<-tryCatch(coloc::coloc.susie(fits[['qtl']],fits[[trait]],p1=1e-4,p2=1e-4,p12=prior),error=function(e)e)
  if(inherits(ans,'error')){comp[[length(comp)+1L]]<-data.frame(trait=trait,p12=prior,status='ERROR',detail=conditionMessage(ans));next}
  if(is.null(ans$summary) || nrow(ans$summary)==0){comp[[length(comp)+1L]]<-data.frame(trait=trait,p12=prior,status='NO_CREDIBLE_SET_PAIR',detail='');next}
  df<-as.data.frame(ans$summary);df$trait<-trait;df$p12<-prior;write.table(df,file.path(p,paste0('coloc_',trait,'_',format(prior,scientific=TRUE),'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  comp[[length(comp)+1L]]<-data.frame(trait=trait,p12=prior,status='EXPLORATORY',detail=paste(nrow(df),'credible-set pairs'))
 }
}
write.table(do.call(rbind,comp),file.path(p,'coloc_status.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
writeLines(c(paste('R:',R.version.string),paste('susieR:',as.character(packageVersion('susieR'))),paste('coloc:',as.character(packageVersion('coloc'))),paste('BLAS:',extSoftVersion()[['BLAS']])),file.path(p,'runtime.txt'))

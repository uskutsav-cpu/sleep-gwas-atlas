.libPaths(c('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
stopifnot(requireNamespace('susieR',quietly=TRUE),requireNamespace('coloc',quietly=TRUE))
parent<-'outputs/shared-genetics-atlas/variant_context/rs62264764_dense_qtl/multisignal_pilot';out<-file.path(parent,'outlier_sensitivity')
v<-read.delim(file.path(parent,'variants.tsv'),stringsAsFactors=FALSE);ns<-read.delim(file.path(parent,'sample_sizes.tsv'),stringsAsFactors=FALSE);k<-read.delim(file.path(parent,'kriging_diagnostics.tsv'),stringsAsFactors=FALSE);R<-as.matrix(read.delim(gzfile(file.path(parent,'ld_matrix.tsv.gz')),header=FALSE));stopifnot(identical(v$variant,ns$variant),nrow(R)==nrow(v),nrow(v)==1548)
flag<-unique(k$variant[abs(k$z_std_diff)>5 & k$fit %in% c('qtl','insomnia','adhd')]);keep<-!v$variant %in% flag;stopifnot(all(flag %in% v$variant),sum(keep)>500)
write.table(k[k$variant %in% flag,],file.path(out,'removed_diagnostics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
v<-v[keep,];ns<-ns[keep,];R<-R[keep,keep];dimnames(R)<-list(v$variant,v$variant)
write.table(v['variant'],file.path(out,'retained_variants.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
base<-list()
for(trait in c('insomnia','adhd')){
 l1<-approx.bf.estimates(v$qtl_beta/v$qtl_se,v$qtl_se^2,'quant')$lABF;l2<-approx.bf.estimates(v[[paste0(trait,'_beta')]]/v[[paste0(trait,'_se')]],v[[paste0(trait,'_se')]]^2,'cc')$lABF
 for(p12 in c(1e-5,1e-6)){z<-combine.abf(l1,l2,1e-4,1e-4,p12,quiet=TRUE);base[[length(base)+1L]]<-data.frame(trait=trait,n=nrow(v),p12=p12,H0=z[1],H1=z[2],H2=z[3],H3=z[4],H4=z[5],row.names=NULL)}
}
write.table(do.call(rbind,base),file.path(out,'single_signal.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
settings<-list(qtl=list(z=v$qtl_beta/v$qtl_se,n=233),insomnia=list(z=v$insomnia_beta/v$insomnia_se,n=313750),adhd=list(z=v$adhd_beta/v$adhd_se,n=128214),adhd_minN=list(z=v$adhd_beta/v$adhd_se,n=64699))
set.seed(20261027);fits<-list();report<-list();credible<-list()
for(name in names(settings)){
 d<-settings[[name]];names(d$z)<-v$variant;notes<-character();cat('starting',name,'\n');flush.console()
 fit<-tryCatch(withCallingHandlers(susieR::susie_rss(z=d$z,R=R,n=d$n,L=10,estimate_residual_variance=FALSE,max_iter=100,tol=1e-3),warning=function(w){notes<<-c(notes,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(fit,'error')){report[[length(report)+1L]]<-data.frame(fit=name,N=d$n,status='ERROR',converged=FALSE,iterations=NA,credible_sets=NA,warnings=paste(c(notes,conditionMessage(fit)),collapse=' | '));next}
 fits[[name]]<-fit;saveRDS(fit,file.path(out,paste0('fit_',name,'.rds')));ncs<-length(fit$sets$cs);report[[length(report)+1L]]<-data.frame(fit=name,N=d$n,status='OK',converged=isTRUE(fit$converged),iterations=fit$niter,credible_sets=ncs,warnings=paste(notes,collapse=' | '));cat('finished',name,'converged',isTRUE(fit$converged),'sets',ncs,'\n');flush.console()
 if(ncs>0)for(i in seq_len(ncs)){ix<-fit$sets$cs[[i]];purity<-fit$sets$purity[i,,drop=FALSE];credible[[length(credible)+1L]]<-data.frame(fit=name,credible_set=names(fit$sets$cs)[i],size=length(ix),min_abs_corr=purity[['min.abs.corr']],mean_abs_corr=purity[['mean.abs.corr']],variants=paste(v$variant[ix],collapse=','))}
}
write.table(do.call(rbind,report),file.path(out,'fit_summary.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
if(length(credible))write.table(do.call(rbind,credible),file.path(out,'credible_sets.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
comparison<-list()
for(t in c('insomnia','adhd','adhd_minN'))for(prior in c(1e-5,1e-6)){
 if(is.null(fits[['qtl']])||is.null(fits[[t]])||!isTRUE(fits[['qtl']]$converged)||!isTRUE(fits[[t]]$converged)){comparison[[length(comparison)+1L]]<-data.frame(trait=t,p12=prior,status='FIT_NOT_ESTIMABLE');next}
 fit<-tryCatch(coloc::coloc.susie(fits[['qtl']],fits[[t]],p1=1e-4,p2=1e-4,p12=prior),error=function(e)e)
 if(inherits(fit,'error')){comparison[[length(comparison)+1L]]<-data.frame(trait=t,p12=prior,status=paste('ERROR',conditionMessage(fit)));next}
 if(is.null(fit$summary)||nrow(fit$summary)==0){comparison[[length(comparison)+1L]]<-data.frame(trait=t,p12=prior,status='NO_CREDIBLE_SET_PAIR');next}
 df<-as.data.frame(fit$summary);df$trait<-t;df$p12<-prior;write.table(df,file.path(out,paste0('coloc_',t,'_',format(prior,scientific=TRUE),'.tsv')),sep='\t',quote=FALSE,row.names=FALSE);comparison[[length(comparison)+1L]]<-data.frame(trait=t,p12=prior,status=paste('EXPLORATORY',nrow(df),'pairs'))
}
write.table(do.call(rbind,comparison),file.path(out,'coloc_status.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
diag<-list()
for(t in c('qtl','insomnia','adhd')){
 d<-settings[[t]];names(d$z)<-v$variant;cat('LD diagnostic',t,'\n');flush.console();lambda<-tryCatch(susieR::estimate_s_rss(d$z,R,n=d$n),error=function(e)NA_real_)
 if(!is.finite(lambda)){diag[[length(diag)+1L]]<-data.frame(trait=t,lambda=NA,max_abs_z_std_diff=NA,count_gt3=NA,count_gt5=NA);next}
 x<-suppressWarnings(susieR::kriging_rss(d$z,R,n=d$n,s=lambda))$conditional_dist;diag[[length(diag)+1L]]<-data.frame(trait=t,lambda=lambda,max_abs_z_std_diff=max(abs(x$z_std_diff)),count_gt3=sum(abs(x$z_std_diff)>3),count_gt5=sum(abs(x$z_std_diff)>5))
}
write.table(do.call(rbind,diag),file.path(out,'kriging_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)

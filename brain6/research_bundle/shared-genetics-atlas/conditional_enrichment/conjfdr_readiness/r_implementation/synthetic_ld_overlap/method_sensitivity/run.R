#!/usr/bin/env Rscript
.libPaths(c('/Volumes/Brain6Work/cfdr-order-empty-lib','/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
suppressPackageStartupMessages(library(cfdr.pleio));suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(Matrix))
out<-'outputs/shared-genetics-atlas/conditional_enrichment/conjfdr_readiness/r_implementation/synthetic_ld_overlap/method_sensitivity'
root<-'/Volumes/Brain6Work/cfdr-synthetic-ld-overlap-20261027';refdir<-file.path(root,'reference');dest<-file.path(root,'method_sensitivity_seed1rho0');stopifnot(!dir.exists(dest))
ref<-readRDS(file.path(refdir,'all_chr_perVariant.rds'));n<-nrow(ref);bs<-10L;nb<-n/bs;block<-rep(seq_len(nb),each=bs)
set.seed(69001);e1<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n);e2<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n);causal<-sample(nb,round(.05*nb));z1<-e1+4*(block %in% causal);z2<-e2
t1<-data.table(SNP=ref$SNP,BETA=z1/100,PVAL=2*pnorm(-abs(z1)));t2<-data.table(SNP=ref$SNP,BETA=z2/100,PVAL=2*pnorm(-abs(z2)))
options(warn=1);o<-cfdr_pleio$new();o$init_data(t1,t2,refdat=refdata_location(refdir),local_refdat_path=dest,verbose=FALSE);o$initialize_pruning_index(n_iter=5,seed=516,verbose=FALSE)
stopifnot(all(colSums(o$index)==nb));for(j in 1:5)stopifnot(all(rowsum(as.integer(o$index[,j]),block)==1L))
configs<-list(default=list(smooth=TRUE,adjust=TRUE,correct_gc=TRUE),no_adjust=list(smooth=TRUE,adjust=FALSE,correct_gc=TRUE),no_smoothing=list(smooth=FALSE,adjust=TRUE,correct_gc=TRUE),no_genomic_correction=list(smooth=TRUE,adjust=TRUE,correct_gc=FALSE))
summary<-list();default_hits<-NULL
for(name in names(configs)){
 cfg<-configs[[name]]
 res<-tryCatch({for(i in 1:2)do.call(o$calculate_cond_fdr,c(list(fdr_trait=i,verbose=FALSE),cfg));o$get_trait_results()},error=function(e)e)
 if(inherits(res,'error')){summary[[name]]<-data.frame(config=name,status='ERROR',message=conditionMessage(res),finite=NA,nonfinite=NA,hits=NA,min=NA);next}
 v<-res$conj_fdr;ok<-is.finite(v)
 summary[[name]]<-data.frame(config=name,status='OK',message='',finite=sum(ok),nonfinite=sum(!ok),hits=sum(v[ok]<=.05),min=if(any(ok))min(v[ok]) else NA)
 if(name=='default'){
  stopifnot(sum(v[ok]<=.05)==1,abs(min(v[ok])-0.0424681369491532)<1e-12)
  ix<-which(ok & v<=.05);m<-match(res$SNP[ix],ref$SNP);stopifnot(all(!is.na(m)))
  default_hits<-data.table(SNP=res$SNP[ix],p_trait1=t1$PVAL[m],p_trait2=t2$PVAL[m],
   block=block[m],trait1_nonnull=block[m] %in% causal,cfdr12=res$cfdr12[ix],cfdr21=res$cfdr21[ix],conj_fdr=v[ix])
 }
}
write.table(rbindlist(summary,fill=TRUE),file.path(out,'sensitivity.tsv'),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
write.table(default_hits,file.path(out,'default_hits.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('METHOD_SENSITIVITY_COMPLETED\n');print(sessionInfo())

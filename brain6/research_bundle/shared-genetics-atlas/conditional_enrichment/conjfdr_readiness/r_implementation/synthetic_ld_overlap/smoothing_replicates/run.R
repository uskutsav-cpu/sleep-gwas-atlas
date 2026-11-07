#!/usr/bin/env Rscript
.libPaths(c('/Volumes/Brain6Work/cfdr-order-empty-lib','/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
suppressPackageStartupMessages(library(cfdr.pleio));suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(Matrix))
out<-'outputs/shared-genetics-atlas/conditional_enrichment/conjfdr_readiness/r_implementation/synthetic_ld_overlap/smoothing_replicates'
root<-'/Volumes/Brain6Work/cfdr-synthetic-ld-overlap-20261027';refdir<-file.path(root,'reference');work<-file.path(root,'smoothing_replicates_retry');stopifnot(!dir.exists(work));dir.create(work)
ref<-readRDS(file.path(refdir,'all_chr_perVariant.rds'));n<-nrow(ref);nb<-n/10L;block<-rep(seq_len(nb),each=10L);stopifnot(n==44000)
orig<-fread(file.path(dirname(out),'runs.tsv'))[rho==0];stopifnot(nrow(orig)==10)
rows<-list();options(warn=1)
for(seed in 69001:69010){
 set.seed(seed);e1<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n);e2<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n);causal<-sample(nb,round(.05*nb));z1<-e1+4*(block %in% causal);z2<-e2
 t1<-data.table(SNP=ref$SNP,BETA=z1/100,PVAL=2*pnorm(-abs(z1)));t2<-data.table(SNP=ref$SNP,BETA=z2/100,PVAL=2*pnorm(-abs(z2)))
 o<-cfdr_pleio$new();dest<-file.path(work,paste0('seed',seed));stopifnot(!dir.exists(dest))
 o$init_data(t1,t2,refdat=refdata_location(refdir),local_refdat_path=dest,verbose=FALSE);o$initialize_pruning_index(n_iter=5,seed=516,verbose=FALSE)
 stopifnot(all(colSums(o$index)==nb));for(j in 1:5)stopifnot(all(rowsum(as.integer(o$index[,j]),block)==1L))
 for(smooth in c(TRUE,FALSE)){
  res<-tryCatch({for(i in 1:2)o$calculate_cond_fdr(fdr_trait=i,smooth=smooth,verbose=FALSE);o$get_trait_results()},error=function(e)e)
  if(inherits(res,'error')){row<-data.frame(seed=seed,smooth=smooth,status='ERROR',message=conditionMessage(res),finite=NA,nonfinite=NA,rejections=NA,min_conj=NA)}else{
   v<-res$conj_fdr;ok<-is.finite(v);row<-data.frame(seed=seed,smooth=smooth,status='OK',message='',finite=sum(ok),nonfinite=sum(!ok),rejections=sum(v[ok]<=.05),min_conj=if(any(ok))min(v[ok]) else NA)
   if(smooth){idx<-which(orig[['seed']]==seed);prior<-orig[idx];stopifnot(nrow(prior)==1,row$rejections==prior$rejections,abs(row$min_conj-prior$min_conj)<1e-12)}
  }
  rows[[length(rows)+1]]<-row;write.table(do.call(rbind,rows),file.path(out,'runs.tsv'),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
  cat('completed',seed,'smooth',smooth,'status',row$status,'rejections',row$rejections,'\n');flush.console()
 }
}
cat('ALL_COMPLETED\n');print(sessionInfo())

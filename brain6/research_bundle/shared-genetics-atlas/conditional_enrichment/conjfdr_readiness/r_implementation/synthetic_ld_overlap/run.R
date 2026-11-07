#!/usr/bin/env Rscript
.libPaths(c('/Volumes/Brain6Work/cfdr-order-empty-lib','/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
suppressPackageStartupMessages(library(cfdr.pleio));suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(Matrix))
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1,!dir.exists(args[1]))
work<-args[1];dir.create(work);refdir<-file.path(work,'reference');dir.create(refdir)
out<-'outputs/shared-genetics-atlas/conditional_enrichment/conjfdr_readiness/r_implementation/synthetic_ld_overlap'
nchr<-22L;perchr<-2000L;bs<-10L;n<-nchr*perchr;nb<-n/bs
ref<-data.table(SNP=sprintf('s%06d',seq_len(n)),CHR=rep(seq_len(nchr),each=perchr),
 BP=rep(seq_len(perchr),nchr),A1='A',A2='C',MAF=.2,INTERGENIC=1L)
saveRDS(ref,file.path(refdir,'all_chr_perVariant.rds'))
for(chr in 1:nchr){ids<-ref[CHR==chr,SNP];block<-rep(seq_len(perchr/bs),each=bs)
 ij<-which(outer(block,block,'==') & row(matrix(0,perchr,perchr))!=col(matrix(0,perchr,perchr)),arr.ind=TRUE)
 mat<-sparseMatrix(i=ij[,1],j=ij[,2],x=1,dims=c(perchr,perchr),dimnames=list(ids,ids))
 stopifnot(inherits(mat,'dgCMatrix'),length(mat@x)==perchr*(bs-1))
 saveRDS(mat,file.path(refdir,sprintf('chr%02d_LDpairs.rds',chr)))}
block<-rep(seq_len(nb),each=bs);rows<-list()
for(seed in 1:10){set.seed(69000+seed)
 e1<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n)
 e2<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n)
 causal<-sample(nb,round(.05*nb));z1<-e1+4*(block %in% causal)
 for(rho in c(0,.3)){z2<-rho*e1+sqrt(1-rho^2)*e2
  t1<-data.table(SNP=ref$SNP,BETA=z1/100,PVAL=2*pnorm(-abs(z1)))
  t2<-data.table(SNP=ref$SNP,BETA=z2/100,PVAL=2*pnorm(-abs(z2)))
  o<-cfdr_pleio$new();dest<-file.path(work,paste0('seed',seed,'rho',rho*10))
  stopifnot(!dir.exists(dest));o$init_data(t1,t2,refdat=refdata_location(refdir),local_refdat_path=dest,verbose=FALSE)
  o$initialize_pruning_index(n_iter=5,seed=516,verbose=FALSE)
  stopifnot(identical(dim(o$index),c(n,5L)),all(colSums(o$index)==nb))
  for(j in 1:5) stopifnot(all(rowsum(as.integer(o$index[,j]),block)==1L))
  o$calculate_cond_fdr(fdr_trait=1,verbose=FALSE);o$calculate_cond_fdr(fdr_trait=2,verbose=FALSE)
  res<-o$get_trait_results();stopifnot(nrow(res)==n,all(is.finite(res$conj_fdr)),all(res$conj_fdr>=0 & res$conj_fdr<=1))
  rows[[length(rows)+1]]<-data.frame(seed=69000+seed,rho=rho,n=n,nonnull_blocks=length(causal),
   pruning_iterations=5,rejections=sum(res$conj_fdr<=.05),min_conj=min(res$conj_fdr),
   gc1=o$gc_correct$corrfac[1],gc2=o$gc_correct$corrfac[2])
  write.table(do.call(rbind,rows),file.path(out,'runs.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  cat('completed',seed,rho,'rejections',rows[[length(rows)]]$rejections,'\n');flush.console()
 }}
print(sessionInfo())

#!/usr/bin/env Rscript
# Independently regenerate seed 69001, rho=0, preserving diagnostic warnings.
.libPaths(c('/Volumes/Brain6Work/cfdr-order-empty-lib','/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
suppressPackageStartupMessages(library(cfdr.pleio));suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(Matrix))
root<-'/Volumes/Brain6Work/cfdr-synthetic-ld-overlap-20261027';refdir<-file.path(root,'reference');dest<-file.path(root,'warning_recheck_seed1rho0');stopifnot(!dir.exists(dest))
ref<-readRDS(file.path(refdir,'all_chr_perVariant.rds'));n<-nrow(ref);bs<-10L;nb<-n/bs;block<-rep(seq_len(nb),each=bs)
set.seed(69001);e1<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n);e2<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n);causal<-sample(nb,round(.05*nb));z1<-e1+4*(block %in% causal);z2<-e2
t1<-data.table(SNP=ref$SNP,BETA=z1/100,PVAL=2*pnorm(-abs(z1)));t2<-data.table(SNP=ref$SNP,BETA=z2/100,PVAL=2*pnorm(-abs(z2)))
options(warn=1);o<-cfdr_pleio$new();o$init_data(t1,t2,refdat=refdata_location(refdir),local_refdat_path=dest,verbose=FALSE);o$initialize_pruning_index(n_iter=5,seed=516,verbose=FALSE)
stopifnot(all(colSums(o$index)==nb));for(j in 1:5)stopifnot(all(rowsum(as.integer(o$index[,j]),block)==1L))
o$calculate_cond_fdr(fdr_trait=1,verbose=FALSE);o$calculate_cond_fdr(fdr_trait=2,verbose=FALSE);r<-o$get_trait_results()
stopifnot(sum(r$conj_fdr<=.05)==1,abs(min(r$conj_fdr)-0.0424681369491532)<1e-12)
cat('RECHECK_PASSED min_conj=',min(r$conj_fdr),' rejections=',sum(r$conj_fdr<=.05),'\n');print(sessionInfo())

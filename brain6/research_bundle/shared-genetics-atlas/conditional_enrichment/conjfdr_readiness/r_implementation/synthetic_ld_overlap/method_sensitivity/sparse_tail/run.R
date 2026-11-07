#!/usr/bin/env Rscript
.libPaths(c('/Volumes/Brain6Work/cfdr-order-empty-lib','/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
suppressPackageStartupMessages(library(cfdr.pleio));suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(Matrix));suppressPackageStartupMessages(library(jsonlite))
out<-'outputs/shared-genetics-atlas/conditional_enrichment/conjfdr_readiness/r_implementation/synthetic_ld_overlap/method_sensitivity/sparse_tail'
root<-'/Volumes/Brain6Work/cfdr-synthetic-ld-overlap-20261027';refdir<-file.path(root,'reference');dest<-file.path(root,'sparse_tail_seed1rho0');stopifnot(!dir.exists(dest))
ref<-readRDS(file.path(refdir,'all_chr_perVariant.rds'));n<-nrow(ref);nb<-n/10L;block<-rep(seq_len(nb),each=10L)
set.seed(69001);e1<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n);e2<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n);causal<-sample(nb,round(.05*nb));z1<-e1+4*(block %in% causal);z2<-e2
t1<-data.table(SNP=ref$SNP,BETA=z1/100,PVAL=2*pnorm(-abs(z1)));t2<-data.table(SNP=ref$SNP,BETA=z2/100,PVAL=2*pnorm(-abs(z2)))
options(warn=1);o<-cfdr_pleio$new();o$init_data(t1,t2,refdat=refdata_location(refdir),local_refdat_path=dest,verbose=FALSE);o$initialize_pruning_index(n_iter=5,seed=516,verbose=FALSE)
stopifnot(all(colSums(o$index)==nb));for(j in 1:5)stopifnot(all(rowsum(as.integer(o$index[,j]),block)==1L))
ix<-which(o$trait_data$SNP=='s021341');stopifnot(length(ix)==1)
fields<-c('LOG10PVAL1','LOG10PVAL2','LOG10PVAL1corr','LOG10PVAL2corr');stopifnot(all(fields %in% names(o$trait_data)))
target<-o$trait_data[ix];counts<-list()
for(j in 1:5){sub<-o$trait_data[o$index[,j]];a<-sub$LOG10PVAL1corr >= target$LOG10PVAL1corr;b<-sub$LOG10PVAL2corr >= target$LOG10PVAL2corr
 counts[[j]]<-data.frame(iteration=j,retained=nrow(sub),trait1_tail=sum(a),trait2_tail=sum(b),both=sum(a & b),target_retained=o$index[ix,j])}
outputs<-list();for(smooth in c(TRUE,FALSE)){
 for(i in 1:2)o$calculate_cond_fdr(fdr_trait=i,smooth=smooth,verbose=FALSE)
 r<-o$get_trait_results();stopifnot(r$SNP[ix]=='s021341');outputs[[as.character(smooth)]]<-list(cfdr12=r$cfdr12[ix],cfdr21=r$cfdr21[ix],conj=r$conj_fdr[ix])}
stopifnot(abs(outputs[['TRUE']]$conj-0.0424681369491532)<1e-12)
summary<-list(snp='s021341',raw_p1=10^(-target$LOG10PVAL1),raw_p2=10^(-target$LOG10PVAL2),
 corrected_p1=10^(-target$LOG10PVAL1corr),corrected_p2=10^(-target$LOG10PVAL2corr),
 outputs=outputs,exact_tail_counts=do.call(rbind,counts),note='post-selected exact threshold counts, not package grid bins')
write_json(summary,file.path(out,'results.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
cat('SPARSE_TAIL_COMPLETED\n');print(sessionInfo())

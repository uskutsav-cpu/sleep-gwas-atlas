#!/usr/bin/env Rscript
.libPaths(c('/Volumes/Brain6Work/cfdr-order-empty-lib','/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
suppressPackageStartupMessages(library(cfdr.pleio));suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(Matrix))
out<-'outputs/shared-genetics-atlas/conditional_enrichment/conjfdr_readiness/r_implementation/synthetic_ld_overlap/mixed_signal'
refdir<-'/Volumes/Brain6Work/cfdr-synthetic-ld-overlap-20261027/reference'
old<-'/Volumes/Extreme SSD/brain6-cfdr-mixed-signal-20261028'
work<-'/Volumes/Brain6Work/cfdr-mixed-resume-20261028';stopifnot(!dir.exists(work));dir.create(work)
cache<-file.path(work,'reference_cache');ref<-readRDS(file.path(refdir,'all_chr_perVariant.rds'))
n<-nrow(ref);nb<-n/10L;stopifnot(n==44000,nb==4400);block<-rep(seq_len(nb),each=10L)
oldtab<-read.delim(file.path(out,'runs.tsv'));stopifnot(nrow(oldtab)==20)
rows<-list();cachehash<-NULL;options(warn=1)
summarize<-function(method,score,truth,seed,rho){
 stopifnot(length(score)==n,all(is.finite(score)),all(score>=0),all(score<=1))
 hit<-score<=.05;true_blocks<-unique(block[truth]);detected<-unique(block[hit]);false_blocks<-setdiff(detected,true_blocks)
 data.frame(seed=seed,rho=rho,method=method,markers=sum(hit),false_markers=sum(hit & !truth),
  fdp_markers=sum(hit & !truth)/max(sum(hit),1),blocks=length(detected),false_blocks=length(false_blocks),
  fdp_blocks=length(false_blocks)/max(length(detected),1),shared_blocks=length(true_blocks),
  shared_blocks_detected=length(intersect(detected,true_blocks)),min_score=min(score))
}
for(seed in 71001:71010){
 set.seed(seed);e1<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n);e2<-sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n)
 chosen<-sample(nb,242);shared<-chosen[1:22];only1<-chosen[23:220];only2<-chosen[221:242]
 truth<-block %in% shared;mean1<-4*(block %in% c(shared,only1));mean2<-4*(block %in% c(shared,only2))
 for(rho in c(0,.3)){
  check_old<-seed==71001 && rho==0
  if(!check_old && any(oldtab$seed==seed & oldtab$rho==rho))next
  z1<-e1+mean1;z2<-rho*e1+sqrt(1-rho^2)*e2+mean2;p1<-2*pnorm(-abs(z1));p2<-2*pnorm(-abs(z2));u<-pmax(p1,p2)
  o<-cfdr_pleio$new();o$init_data(data.table(SNP=ref$SNP,BETA=z1/100,PVAL=p1),data.table(SNP=ref$SNP,BETA=z2/100,PVAL=p2),
   refdat=refdata_location(refdir),local_refdat_path=cache,verbose=FALSE)
  paths<-sort(list.files(cache,pattern='\\.rds$',full.names=TRUE));stopifnot(length(paths)==23)
  hashes<-unname(tools::md5sum(paths));if(is.null(cachehash))cachehash<-hashes else stopifnot(identical(cachehash,hashes))
  o$initialize_pruning_index(n_iter=5,seed=516,verbose=FALSE)
  stopifnot(all(colSums(o$index)==nb));for(j in 1:5)stopifnot(all(rowsum(as.integer(o$index[,j]),block)==1L))
  scores<-list();for(smooth in c(TRUE,FALSE)){
   for(i in 1:2)o$calculate_cond_fdr(fdr_trait=i,smooth=smooth,verbose=FALSE)
   res<-o$get_trait_results();stopifnot(nrow(res)==n,identical(as.character(res$SNP),as.character(ref$SNP)))
   scores[[if(smooth)'smooth' else 'unsmoothed']]<-res$conj_fdr
  }
  scores$maxp_bonferroni<-p.adjust(u,'bonferroni');scores$maxp_BY<-p.adjust(u,'BY')
  if(check_old){
   prior<-readRDS(file.path(old,'scores_71001_rho0.rds'))
   stopifnot(identical(z1,prior$z1),identical(z2,prior$z2),identical(truth,prior$truth))
   for(key in names(scores))stopifnot(isTRUE(all.equal(scores[[key]],prior$scores[[key]],tolerance=1e-12)))
   cat('PASS original seed71001 rho0 reproduced; cache hashes recorded\n');flush.console();next
  }
  for(key in names(scores))rows[[length(rows)+1]]<-summarize(key,scores[[key]],truth,seed,rho)
  saveRDS(list(seed=seed,rho=rho,z1=z1,z2=z2,block=block,truth=truth,shared=shared,only1=only1,only2=only2,scores=scores),
   file.path(work,paste0('scores_',seed,'_rho',rho*10,'.rds')))
  write.table(do.call(rbind,rows),file.path(out,'resume_runs.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  cat('completed seed',seed,'rho',rho,'\n');flush.console()
 }
}
write.table(data.frame(path=basename(paths),md5=cachehash),file.path(out,'reference_cache_md5.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('ALL_COMPLETED\n');print(sessionInfo())

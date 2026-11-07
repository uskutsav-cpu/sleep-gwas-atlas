#!/usr/bin/env Rscript
root<-'outputs/shared-genetics-atlas/conditional_enrichment/conjfdr_readiness/r_implementation/synthetic_ld_overlap/mixed_signal'
old<-'/Volumes/Extreme SSD/brain6-cfdr-mixed-signal-20261028';new<-'/Volumes/Brain6Work/cfdr-mixed-resume-20261028'
tab<-rbind(read.delim(file.path(root,'runs.tsv')),read.delim(file.path(root,'resume_runs.tsv')))
stopifnot(nrow(tab)==80,!anyDuplicated(tab[c('seed','rho','method')]))
files<-c(list.files(old,pattern='^scores_.*\\.rds$',full.names=TRUE),list.files(new,pattern='^scores_.*\\.rds$',full.names=TRUE));stopifnot(length(files)==20)
for(path in files){
 x<-readRDS(path);n<-length(x$z1);stopifnot(n==44000,length(x$shared)==22,length(x$only1)==198,length(x$only2)==22,
  length(unique(c(x$shared,x$only1,x$only2)))==242,identical(x$truth,x$block %in% x$shared))
 pmax2<-pmax(2*pnorm(-abs(x$z1)),2*pnorm(-abs(x$z2)))
 stopifnot(isTRUE(all.equal(x$scores$maxp_BY,p.adjust(pmax2,'BY'),tolerance=1e-13)),
  isTRUE(all.equal(x$scores$maxp_bonferroni,p.adjust(pmax2,'bonferroni'),tolerance=1e-13)))
 for(method in names(x$scores)){
  score<-x$scores[[method]];hit<-score<=.05;row<-tab[tab$seed==x$seed & tab$rho==x$rho & tab$method==method,];stopifnot(nrow(row)==1)
  db<-unique(x$block[hit]);false<-setdiff(db,x$shared)
  stopifnot(row$markers==sum(hit),row$false_markers==sum(hit & !x$truth),row$blocks==length(db),
   row$false_blocks==length(false),row$shared_blocks_detected==length(intersect(db,x$shared)),
   abs(row$fdp_markers-sum(hit & !x$truth)/max(1,sum(hit)))<1e-12,
   abs(row$fdp_blocks-length(false)/max(1,length(db)))<1e-12)
 }
}
summary<-do.call(rbind,lapply(split(tab,list(tab$rho,tab$method),drop=TRUE),function(x)
 data.frame(rho=x$rho[1],method=x$method[1],replicates=nrow(x),mean_marker_fdp=mean(x$fdp_markers),
  mcse_marker_fdp=sd(x$fdp_markers)/sqrt(nrow(x)),mean_block_fdp=mean(x$fdp_blocks),
  mean_shared_blocks_detected=mean(x$shared_blocks_detected),runs_with_false=sum(x$false_markers>0),
  mean_markers=mean(x$markers))))
write.table(summary,file.path(root,'complete_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS: 20 stored datasets, 80 rows, truth labels, marker/block counts and max-p adjustments verified.\n');print(summary,row.names=FALSE)

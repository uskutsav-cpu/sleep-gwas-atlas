base <- '/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support'
ref <- '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz'
out <- 'outputs/shared-genetics-atlas/conditional_enrichment/block_sensitivity'
s <- read.delim(gzfile(file.path(base,'S.sumstats.gz')),stringsAsFactors=FALSE)
f <- read.delim(gzfile(file.path(base,'FI.sumstats.gz')),stringsAsFactors=FALSE)
t <- read.delim(gzfile(file.path(base,'FFS.sumstats.gz')),stringsAsFactors=FALSE)
r <- read.delim(gzfile(ref),stringsAsFactors=FALSE)
stopifnot(nrow(s)==1138262,!anyDuplicated(s$SNP),identical(s$SNP,f$SNP),identical(s$SNP,t$SNP),identical(s$A1,f$A1),identical(s$A2,f$A2),identical(s$A1,t$A1),identical(s$A2,t$A2))
dup<-duplicated(r$SNP)|duplicated(r$SNP,fromLast=TRUE);j<-match(s$SNP,r$SNP)
valid<-!is.na(j);valid[valid]<-!dup[j[valid]] & ((s$A1[valid]==r$A1[j[valid]] & s$A2[valid]==r$A2[j[valid]]) | (s$A1[valid]==r$A2[j[valid]] & s$A2[valid]==r$A1[j[valid]])) & r$CHR[j[valid]] %in% 1:22 & !is.na(r$BP[j[valid]]) & r$BP[j[valid]]>0
cat('support',nrow(s),'matched',sum(valid),'held',sum(!valid),'\n');stopifnot(sum(valid)>1e6)
chr<-r$CHR[j[valid]];bp<-r$BP[j[valid]];zs<-s$Z[valid];zf<-f$Z[valid];zt<-t$Z[valid];ps<-2*pnorm(-abs(zs));block_rows<-list(); summary_rows<-list();count<-0L;k<-0L
for(width in c(5000000,10000000)) {
 block<-paste(chr,floor((bp-1)/width),sep=':');group<-match(block,unique(block));levels<-unique(block)
 for(pair in c('FI','FFS')) for(direction in c('insomnia_given_frailty','frailty_given_insomnia')) {
  pf<-2*pnorm(-abs(if(pair=='FI')zf else zt));target<-if(direction=='insomnia_given_frailty')ps else pf;condition<-if(direction=='insomnia_given_frailty')pf else ps
  is_target<-target<=.01;is_cond<-condition<=.001
  n<-tabulate(group,nbins=length(levels));nt<-tabulate(group[is_target],nbins=length(levels));nc<-tabulate(group[is_cond],nbins=length(levels));nj<-tabulate(group[is_target & is_cond],nbins=length(levels))
  N<-sum(n);T<-sum(nt);C<-sum(nc);J<-sum(nj);stopifnot(C>0,all(N-n>0),all(C-nc>0))
  baseline<-T/N;conditional<-J/C;delta<-conditional-baseline;leave<-(J-nj)/(C-nc)-(T-nt)/(N-n);B<-length(leave);se<-sqrt((B-1)/B*sum((leave-mean(leave))^2))
  count<-count+1L;summary_rows[[count]]<-data.frame(pair=pair,direction=direction,block_bp=width,n_matched=N,n_blocks=B,baseline=baseline,conditional=conditional,difference=delta,ratio=conditional/baseline,jackknife_se=se,normal_95_lower=delta-1.96*se,normal_95_upper=delta+1.96*se,max_block_fraction_of_condition=max(nc)/C)
  k<-k+1L;block_rows[[k]]<-data.frame(pair=pair,direction=direction,block_bp=width,block=levels,n=n,n_target=nt,n_condition=nc,n_joint=nj)
 }
}
write.table(do.call(rbind,summary_rows),file.path(out,'summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,block_rows),file.path(out,'block_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(do.call(rbind,summary_rows))

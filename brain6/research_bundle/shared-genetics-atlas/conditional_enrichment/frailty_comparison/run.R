base<-'/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation'
ref<-'/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz'
out<-'outputs/shared-genetics-atlas/conditional_enrichment/frailty_comparison'
s<-read.delim(gzfile(file.path(base,'common_support/S.sumstats.gz')),stringsAsFactors=FALSE)
f<-read.delim(gzfile(file.path(base,'common_support/FI.sumstats.gz')),stringsAsFactors=FALSE)
t<-read.delim(gzfile(file.path(base,'common_support/FFS.sumstats.gz')),stringsAsFactors=FALSE)
h<-read.delim(file.path(base,'filtered/HFRS_HM3.tsv'),stringsAsFactors=FALSE)
r<-read.delim(gzfile(ref),stringsAsFactors=FALSE)
stopifnot(nrow(s)==1138262L,!anyDuplicated(s$SNP),identical(s$SNP,f$SNP),identical(s$SNP,t$SNP),identical(s$A1,f$A1),identical(s$A2,f$A2),identical(s$A1,t$A1),identical(s$A2,t$A2))
hrdup<-duplicated(h$SNP)|duplicated(h$SNP,fromLast=TRUE); rrdup<-duplicated(r$SNP)|duplicated(r$SNP,fromLast=TRUE);jh<-match(s$SNP,h$SNP);jr<-match(s$SNP,r$SNP)
valid<-!is.na(jh)&!is.na(jr); idx<-which(valid)
valid[idx]<-!hrdup[jh[idx]] & !rrdup[jr[idx]] & ((s$A1[idx]==h$A1[jh[idx]] & s$A2[idx]==h$A2[jh[idx]]) | (s$A1[idx]==h$A2[jh[idx]] & s$A2[idx]==h$A1[jh[idx]])) & ((s$A1[idx]==r$A1[jr[idx]] & s$A2[idx]==r$A2[jr[idx]]) | (s$A1[idx]==r$A2[jr[idx]] & s$A2[idx]==r$A1[jr[idx]])) & h$CHR[jh[idx]]==r$CHR[jr[idx]] & r$CHR[jr[idx]] %in% 1:22 & !is.na(r$BP[jr[idx]]) & r$BP[jr[idx]]>0 & is.finite(h$P[jh[idx]]) & h$P[jh[idx]]>0 & h$P[jh[idx]]<=1
stopifnot(all(is.finite(s$Z)),all(is.finite(f$Z)),all(is.finite(t$Z)),sum(valid)>1e6)
cat('original',nrow(s),'joint_support',sum(valid),'held',sum(!valid),'\n')
pos<-r$BP[jr[valid]];chr<-r$CHR[jr[valid]]; ps<-2*pnorm(-abs(s$Z[valid]));cond<-list(FI=2*pnorm(-abs(f$Z[valid]))<=.001, FFS=2*pnorm(-abs(t$Z[valid]))<=.001,HFRS=h$P[jh[valid]]<=.001); target<-ps<=.01
summary<-list();blocks<-list();i<-0L;k<-0L
for(width in c(5000000,10000000)){
 block<-paste(chr,floor((pos-1)/width),sep=':');levels<-unique(block);g<-match(block,levels);B<-length(levels);nn<-tabulate(g,nbins=B);na<-tabulate(g[target],nbins=B)
 cc<-lapply(cond,function(v)tabulate(g[v],nbins=B)); jj<-lapply(cond,function(v)tabulate(g[v&target],nbins=B));basep<-sum(na)/sum(nn)
 b<-data.frame(block_bp=width,block=levels,n=nn,n_insomnia=na,FI_cond=cc$FI,FI_joint=jj$FI,FFS_cond=cc$FFS,FFS_joint=jj$FFS,HFRS_cond=cc$HFRS,HFRS_joint=jj$HFRS);k<-k+1L;blocks[[k]]<-b
 for(pair in list(c('FI','FFS'),c('FI','HFRS'),c('FFS','HFRS'))){
  a<-pair[1];b<-pair[2];ca<-sum(cc[[a]]);cb<-sum(cc[[b]]);ja<-sum(jj[[a]]);jb<-sum(jj[[b]]);stopifnot(all(ca-cc[[a]]>0),all(cb-cc[[b]]>0))
  fa<-ja/ca;fb<-jb/cb;delta<-fa-fb;leave<-(ja-jj[[a]])/(ca-cc[[a]])-(jb-jj[[b]])/(cb-cc[[b]]);se<-sqrt((B-1)/B*sum((leave-mean(leave))^2));i<-i+1L
  summary[[i]]<-data.frame(contrast=paste(a,b,sep=' minus '),block_bp=width,matched_n=sum(nn),blocks=B,baseline_insomnia=basep,selected_a=ca,joint_a=ja,fraction_a=fa,selected_b=cb,joint_b=jb,fraction_b=fb,difference=delta,jackknife_se=se,normal_95_lower=delta-1.96*se,normal_95_upper=delta+1.96*se)
 }
}
write.table(do.call(rbind,summary),file.path(out,'summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,blocks),file.path(out,'block_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(do.call(rbind,summary))

args <- commandArgs(trailingOnly=TRUE); if(length(args)!=3)stop('expected pair input directory and output directory')
pair<-args[1]; input<-args[2]; output<-args[3]; if(!pair %in% c('FI','FFS')) stop('invalid trait')
source('/Volumes/Extreme SSD/brain6-work/PLACO-v0.2.0/PLACO_v0.2.0.R')
s<-read.delim(gzfile(file.path(input,'S.sumstats.gz')),stringsAsFactors=FALSE)
t<-read.delim(gzfile(file.path(input,paste0(pair,'.sumstats.gz'))),stringsAsFactors=FALSE)
stopifnot(nrow(s)==1138262L,nrow(t)==nrow(s),!anyDuplicated(s$SNP),identical(s$SNP,t$SNP),identical(s$A1,t$A1),identical(s$A2,t$A2),all(is.finite(s$Z)),all(is.finite(t$Z)))
settings<-list(FI=c(1.271851342700776,1.4045250827610078,0.08081972497537657),FFS=c(1.2714481467007943,1.5006275551853958,0.052417344112559244));par<-settings[[pair]]
partial<-file.path(output,paste0(pair,'_screen.partial.tsv.gz')); final<-file.path(output,paste0(pair,'_screen.tsv.gz'));stopifnot(!file.exists(partial),!file.exists(final));o<-gzfile(partial,'wt');n<-nrow(s);fail<-0L;zero<-0L;excluded<-0L;found<-0L;elapsed<-proc.time()[3]
for(start in seq.int(1L,n,by=5000L)){
 ix<-start:min(start+4999L,n); z1<-s$Z[ix]; z2<-t$Z[ix]; stat<-vapply(seq_along(ix),function(i){
  if(z1[i]^2>80 || z2[i]^2>80)return(NA_real_)
  if(z1[i]==0 || z2[i]==0)return(1)
  tryCatch({p<-placo.plus(c(z1[i],z2[i]),VarZ=par[1:2],CorZ=par[3],AbsTol=1e-13)$p.placo.plus; if(length(p)==1L && is.finite(p) && p>0 && p<=1)p else NA_real_},error=function(e)NA_real_)
 },numeric(1)); x<-z1^2>80|z2^2>80; status<-ifelse(x,'EXCLUDED_EXTREME_Z',ifelse(is.na(stat),'NUMERICAL_FAILURE','TESTED')); fail<-fail+sum(is.na(stat)&!x);excluded<-excluded+sum(x);zero<-zero+sum(z1==0|z2==0);found<-found+sum(stat<=(0.05/(2*1138262)),na.rm=TRUE)
 dd<-data.frame(SNP=s$SNP[ix],A1=s$A1[ix],A2=s$A2[ix],Z_insomnia=z1,Z_frailty=z2,P_insomnia=2*pnorm(-abs(z1)),P_frailty=2*pnorm(-abs(z2)),P_PLACO_plus=stat,status=status)
 write.table(dd,o,sep='\t',quote=FALSE,row.names=FALSE,col.names=start==1L,na='NA'); flush(o)
 if(start==1L||start%%50000L==1L||max(ix)==n) {cat(pair,'completed',max(ix),'/',n,'failures',fail,'exploratory_threshold_hits',found,'elapsed_s',round(proc.time()[3]-elapsed),'\n');flush.console()}
}
close(o); if(!file.rename(partial,final))stop('result rename failed');cat(pair,'COMPLETE rows=',n,' failures=',fail,' excluded=',excluded,' zeros=',zero,' threshold_hits=',found,' seconds=',round(proc.time()[3]-elapsed),'\n')

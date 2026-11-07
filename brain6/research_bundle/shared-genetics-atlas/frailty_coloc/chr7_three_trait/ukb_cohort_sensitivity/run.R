.libPaths(c('/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()));library(hyprcoloc)
source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p<-'outputs/shared-genetics-atlas/frailty_coloc/chr7_three_trait/ukb_cohort_sensitivity';x<-read.delim('/Volumes/Brain6Work/chr7-three-trait-ukb-sensitivity.tsv');stopifnot(nrow(x)>500,!anyDuplicated(x$rsid))
settings<-data.frame(name=c('default','prior_001','prior_0005','strict_threshold'),prior_c=c(.02,.01,.005,.02),threshold=c(.5,.5,.5,.9));cluster<-list();posterior<-list();models<-list()
for (source in c('meta','UKB')) {
 f<-if(source=='meta') x$fi_beta else x$ukb_fi_beta; fs<-if(source=='meta') x$fi_se else x$ukb_fi_se
 b<-cbind(FI=f,insomnia=x$insomnia_beta,ADHD=x$adhd_beta);s<-cbind(FI=fs,insomnia=x$insomnia_se,ADHD=x$adhd_se)
 stopifnot(all(is.finite(b)),all(is.finite(s)),all(s>0))
 for(k in seq_len(nrow(settings))) {v<-settings[k,];m<-hyprcoloc(b,s,binary.outcomes=c(0,1,1),trait.names=colnames(b),snp.id=x$rsid,prior.1=1e-4,prior.c=v$prior_c,reg.thresh=v$threshold,align.thresh=v$threshold,ind.traits=TRUE,bb.alg=TRUE,snpscores=TRUE);d<-m$results;d$setting<-v$name;d$source<-source;cluster[[length(cluster)+1]]<-d;models[[paste(source,v$name,sep='_')]]<-m}
 score<-function(z,se,type) approx.bf.estimates(z/se,se^2,type)$lABF
 ff<-score(f,fs,'quant');ii<-score(x$insomnia_beta,x$insomnia_se,'cc');aa<-score(x$adhd_beta,x$adhd_se,'cc')
 for(pair in c('insomnia','ADHD')) for(pr in c(1e-5,1e-6)) {v<-combine.abf(ff,if(pair=='insomnia') ii else aa,1e-4,1e-4,pr,quiet=TRUE);stopifnot(abs(sum(v)-1)<1e-12);posterior[[length(posterior)+1]]<-data.frame(source=source,pair=pair,p12=pr,H0=v[1],H1=v[2],H2=v[3],H3=v[4],H4=v[5],row.names=NULL)}
}
write.table(do.call(rbind,cluster),file.path(p,'clusters.tsv'),sep='\t',quote=FALSE,row.names=FALSE);write.table(do.call(rbind,posterior),file.path(p,'pairwise.tsv'),sep='\t',quote=FALSE,row.names=FALSE);saveRDS(models,file.path(p,'models.rds'));writeLines(capture.output(sessionInfo()),file.path(p,'session_info.txt'))

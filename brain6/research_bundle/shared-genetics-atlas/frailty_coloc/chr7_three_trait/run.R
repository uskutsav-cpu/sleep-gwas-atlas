.libPaths(c('/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()));library(hyprcoloc)
source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p<-'outputs/shared-genetics-atlas/frailty_coloc/chr7_three_trait';x<-read.delim('/Volumes/Brain6Work/chr7-three-trait-input.tsv',stringsAsFactors=FALSE)
stopifnot(nrow(x)>=500,!anyDuplicated(x$rsid))
smoke<-hyprcoloc(hyprcoloc::test.betas,hyprcoloc::test.ses,trait.names=paste0('T',1:ncol(hyprcoloc::test.betas)),snp.id=rownames(hyprcoloc::test.betas),ind.traits=TRUE);stopifnot(nrow(smoke$results)==3,identical(as.character(smoke$results$candidate_snp),c('rs11591147','rs12117612','rs7524677')));saveRDS(smoke,file.path(p,'smoke_test.rds'))
b<-cbind(FI=x$fi_beta,insomnia=x$insomnia_beta,ADHD=x$adhd_beta);s<-cbind(FI=x$fi_se,insomnia=x$insomnia_se,ADHD=x$adhd_se);stopifnot(all(is.finite(b)),all(is.finite(s)),all(s>0))
settings<-data.frame(name=c('default','prior_001','prior_0005','strict_threshold'),prior_c=c(.02,.01,.005,.02),threshold=c(.5,.5,.5,.9));res<-list();clusters<-list()
for(k in seq_len(nrow(settings))){v<-settings[k,];cat('\nSETTING:',v$name,'prior.c',v$prior_c,'threshold',v$threshold,'\n');r<-hyprcoloc(b,s,binary.outcomes=c(0,1,1),trait.names=colnames(b),snp.id=x$rsid,prior.1=1e-4,prior.c=v$prior_c,reg.thresh=v$threshold,align.thresh=v$threshold,ind.traits=TRUE,bb.alg=TRUE,snpscores=TRUE);print(r);res[[v$name]]<-r;d<-r$results;d$setting<-v$name;clusters[[v$name]]<-d}
saveRDS(res,file.path(p,'pilot_results.rds'));write.table(settings,file.path(p,'settings.tsv'),sep='\t',row.names=FALSE,quote=FALSE);write.table(do.call(rbind,clusters),file.path(p,'clusters.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
score<-function(z,se,type) approx.bf.estimates(z/se,se^2,type)$lABF
f<-score(x$fi_beta,x$fi_se,'quant');i<-score(x$insomnia_beta,x$insomnia_se,'cc');a<-score(x$adhd_beta,x$adhd_se,'cc');pairs<-list(FI_insomnia=list(f,i),FI_ADHD=list(f,a),insomnia_ADHD=list(i,a));rows<-list();j<-0
for(n in names(pairs))for(prior in c(1e-5,1e-6)){z<-combine.abf(pairs[[n]][[1]],pairs[[n]][[2]],1e-4,1e-4,prior,quiet=TRUE);stopifnot(abs(sum(z)-1)<1e-12);j<-j+1;rows[[j]]<-data.frame(pair=n,n=nrow(x),p12=prior,H0=z[1],H1=z[2],H2=z[3],H3=z[4],H4=z[5],row.names=NULL)}
write.table(do.call(rbind,rows),file.path(p,'pairwise.tsv'),sep='\t',row.names=FALSE,quote=FALSE);print(do.call(rbind,rows));writeLines(capture.output(sessionInfo()),file.path(p,'session_info.txt'))

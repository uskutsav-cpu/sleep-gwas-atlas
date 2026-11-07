.libPaths(c('/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()));library(hyprcoloc)
p <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_dense_qtl/hyprcoloc_pilot'
src <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_dense_qtl/multisignal_pilot/variants.tsv'
cat('R:',R.version.string,'\nHyPrColoc:',as.character(packageVersion('hyprcoloc')),'\n')
smoke<-hyprcoloc(hyprcoloc::test.betas,hyprcoloc::test.ses,trait.names=paste0('T',1:ncol(hyprcoloc::test.betas)),snp.id=rownames(hyprcoloc::test.betas),ind.traits=TRUE)
stopifnot(nrow(smoke$results)==3,identical(as.character(smoke$results$candidate_snp),c('rs11591147','rs12117612','rs7524677')))
saveRDS(smoke,file.path(p,'smoke_test.rds'));cat('Documented example:');print(smoke)
x<-read.delim(src,stringsAsFactors=FALSE);stopifnot(nrow(x)==1548,!anyDuplicated(x$variant))
b<-cbind(expression=x$qtl_beta,insomnia=x$insomnia_beta,adhd=x$adhd_beta)
s<-cbind(expression=x$qtl_se,insomnia=x$insomnia_se,adhd=x$adhd_se)
stopifnot(all(is.finite(b)),all(is.finite(s)),all(s>0))
settings<-data.frame(name=c('default','prior_001','prior_0005','strict_threshold'),prior_c=c(.02,.01,.005,.02),threshold=c(.5,.5,.5,.9));results<-list();rows<-list()
for(k in seq_len(nrow(settings))){v<-settings[k,];cat('\nSETTING:',v$name,'prior.c',v$prior_c,'threshold',v$threshold,'\n')
 result<-hyprcoloc(b,s,binary.outcomes=c(0,1,1),trait.names=colnames(b),snp.id=x$variant,prior.1=1e-4,prior.c=v$prior_c,reg.thresh=v$threshold,align.thresh=v$threshold,ind.traits=TRUE,bb.alg=TRUE,snpscores=TRUE)
 print(result);results[[v$name]]<-result;d<-result$results;d$setting<-v$name;rows[[v$name]]<-d}
saveRDS(results,file.path(p,'pilot_results.rds'));write.table(settings,file.path(p,'settings.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,rows),file.path(p,'clusters.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(p,'session_info.txt'))

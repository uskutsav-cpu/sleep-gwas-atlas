.libPaths(c('/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()));library(hyprcoloc);source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p<-'outputs/shared-genetics-atlas/frailty_coloc/chr7_three_trait/precision_sensitivity';x<-read.delim('/Volumes/Brain6Work/chr7-three-trait-input.tsv',stringsAsFactors=FALSE);stopifnot(nrow(x)==1055)
set.seed(20261027);res<-list();compute<-function(bfi,sfi,n){stopifnot(all(sfi>0));b<-cbind(FI=bfi,insomnia=x$insomnia_beta,ADHD=x$adhd_beta);s<-cbind(FI=sfi,insomnia=x$insomnia_se,ADHD=x$adhd_se)
 h<-hyprcoloc(b,s,binary.outcomes=c(0,1,1),trait.names=colnames(b),snp.id=x$rsid,prior.1=1e-4,prior.c=.02,reg.thresh=.5,align.thresh=.5,ind.traits=TRUE,bb.alg=TRUE,snpscores=FALSE)$results
 sig<-function(beta,se,type)approx.bf.estimates(beta/se,se^2,type)$lABF
 f<-sig(bfi,sfi,'quant');i<-sig(x$insomnia_beta,x$insomnia_se,'cc');a<-sig(x$adhd_beta,x$adhd_se,'cc')
 h4<-function(z,p12)unname(combine.abf(f,z,1e-4,1e-4,p12,quiet=TRUE)[5])
 data.frame(replicate=n,cluster=paste(h$traits,collapse=' | '),posterior=paste(h$posterior_prob,collapse=' | '),candidate=paste(h$candidate_snp,collapse=' | '),FI_insomnia_1e5=h4(i,1e-5),FI_insomnia_1e6=h4(i,1e-6),FI_ADHD_1e5=h4(a,1e-5),FI_ADHD_1e6=h4(a,1e-6))}
res[[1]]<-compute(x$fi_beta,x$fi_se,0)
for(j in 1:200){res[[j+1]]<-compute(x$fi_beta+runif(nrow(x),-.00005,.00005),x$fi_se+runif(nrow(x),-.00005,.00005),j);if(j%%50==0)cat('completed',j,'\n')}
all<-do.call(rbind,res);write.table(all,file.path(p,'replicates.tsv'),sep='\t',quote=FALSE,row.names=FALSE);cat('Cluster frequencies:\n');print(table(all$cluster[-1]));cat('Candidate frequencies:\n');print(sort(table(all$candidate[-1]),decreasing=TRUE));print(all[1,]);for(v in c('FI_insomnia_1e5','FI_insomnia_1e6','FI_ADHD_1e5','FI_ADHD_1e6'))cat(v,'range',range(all[[v]][-1]),'\n');cat('Cluster posterior range',range(as.numeric(all$posterior[-1]),na.rm=TRUE),'\n')

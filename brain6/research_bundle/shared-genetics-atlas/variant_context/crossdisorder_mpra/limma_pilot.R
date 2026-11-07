# Direct matrix implementation of the version 1.24.0 mpra::mpralm corr_groups
# steps for unique element IDs, using limma. Provisional: publication QC and
# independent comparison against the published results table remain unresolved.
suppressPackageStartupMessages(library(limma))
root <- '/Volumes/Extreme SSD/brain6-crossdisorder-mpra'
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args)) args[1] else 'v124'
filter <- if(length(args)>1) args[2] else 'none'
stopifnot(mode %in% c('v124','v134'), filter %in% c('none','all_listed','prop30'))
dna <- read.csv(gzfile(file.path(root,'dnaMatrix.csv.gz')),check.names=FALSE)
rna <- read.csv(gzfile(file.path(root,'rnaMatrix.csv.gz')),check.names=FALSE)
stopifnot(identical(dna$variant,rna$variant),!anyDuplicated(dna$variant),ncol(dna)==45,ncol(rna)==45)
if(filter!='none') { ex<-read.csv(gzfile(file.path(root,'excludelist.csv.gz'))); if(filter=='prop30') ex<-ex[ex$prop>.3,]; exids<-sub('_',':',ex$variant,fixed=TRUE); keep<-!dna$variant %in% exids; dna<-dna[keep,];rna<-rna[keep,] }
ord<-order(dna$variant);ids<-dna$variant[ord];d<-as.matrix(dna[ord,-1]);r<-as.matrix(rna[ord,-1]);storage.mode(d)<-'numeric';storage.mode(r)<-'numeric'
rownames(d)<-rownames(r)<-ids;colnames(d)<-colnames(r)<-c(paste0('Replicate ',1:22,' REF'),paste0('Replicate ',1:22,' ALT'))
block<-rep(1:22,2); batch<-factor(rep(c(rep(1,10),rep(2,12)),2)); alt<-rep(c(FALSE,TRUE),each=22);design<-model.matrix(~alt+batch)
# mpra 1.24.0 .fit_corr calls normalize_counts(object) WITHOUT a block.
ds<-colSums(d,na.rm=TRUE); rs<-colSums(r,na.rm=TRUE)
if(mode=='v134') { ds<-tapply(ds,block,sum,na.rm=TRUE)[block]; rs<-tapply(rs,block,sum,na.rm=TRUE)[block] }
d<-round(sweep(d,2,ds,'/')*1e7);r<-round(sweep(r,2,rs,'/')*1e7)
logdna<-log2(d+1); logr<-log2(r+1)-logdna
unweighted<-lmFit(logr,design);lo<-lowess(rowMeans(logdna,na.rm=TRUE),sqrt(unweighted$sigma),f=.4);fun<-approxfun(lo,rule=2)
weights<-matrix(1/fun(logdna)^4,nrow=nrow(logdna),dimnames=dimnames(logdna))
corfit<-duplicateCorrelation(logr,design=design,ndups=1,block=block)
elist<-new('EList',list(E=logr,weights=weights,design=design))
fit<-eBayes(lmFit(elist,design,block=block,correlation=corfit$consensus))
tab<-topTable(fit,coef=2,number=Inf);tab$variant<-rownames(tab)
write.table(tab,file.path(root,paste0('limma_pilot_',mode,if(filter=='none') '' else paste0('_',filter),'_all.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(paste('consensus correlation',corfit$consensus),paste('rows',nrow(tab)),capture.output(sessionInfo())),file.path(root,paste0('limma_pilot_',mode,if(filter=='none') '' else paste0('_',filter),'_session.txt')))

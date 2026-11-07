# Adapted from Lee et al. public crossdisorder-MPRA pipeline at commit
# 0d6a07b5925031a543b49733d2f9aed17ee69c42. IDs are read from the actual
# variant column; no exclusion or active-element filter is applied here.
suppressPackageStartupMessages(library(mpra))
root <- '/Volumes/Extreme SSD/brain6-crossdisorder-mpra'
dna <- read.csv(gzfile(file.path(root,'dnaMatrix.csv.gz')),check.names=FALSE)
rna <- read.csv(gzfile(file.path(root,'rnaMatrix.csv.gz')),check.names=FALSE)
stopifnot(identical(dna$variant,rna$variant),!anyDuplicated(dna$variant),ncol(dna)==45,ncol(rna)==45)
ids <- dna$variant
d <- as.matrix(dna[,-1]);r <- as.matrix(rna[,-1]);storage.mode(d)<-'numeric';storage.mode(r)<-'numeric'
colnames(d)<-colnames(r)<-c(paste0('Replicate ',1:22,' REF'),paste0('Replicate ',1:22,' ALT'))
batch<-factor(rep(c(rep(1,10),rep(2,12)),2));alt<-rep(c(FALSE,TRUE),each=22);design<-model.matrix(~alt+batch)
set<-MPRASet(DNA=d,RNA=r,eid=ids,eseq=NULL,barcode=NULL)
fit<-mpralm(object=set,design=design,plot=FALSE,aggregate='none',normalize=TRUE,block=rep(1:22,2),model_type='corr_groups')
tab<-topTable(fit,coef=2,number=Inf);tab$variant<-rownames(tab)
write.table(tab,file.path(root,'provisional_allelic_all.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),paste('Rows:',nrow(tab))),file.path(root,'reproduction_session.txt'))

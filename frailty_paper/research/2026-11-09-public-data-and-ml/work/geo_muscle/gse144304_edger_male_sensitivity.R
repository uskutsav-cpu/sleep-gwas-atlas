.libPaths(c('work/rna_Rlib',.libPaths()));suppressPackageStartupMessages(library(edgeR))
meta<-readLines(gzfile('work/geo_muscle/GSE144304.txt.gz'))
parse_row<-function(line) as.character(read.delim(text=line,header=FALSE,quote='"',check.names=FALSE,stringsAsFactors=FALSE)[1,-1])
title<-parse_row(meta[grep('^!Sample_title\t',meta)])
chars<-lapply(meta[grep('^!Sample_characteristics_ch1\t',meta)],parse_row)
extract<-function(key){x<-Filter(function(v) all(startsWith(v,paste0(key,': '))),chars);stopifnot(length(x)==1);sub(paste0('^',key,': '),'',x[[1]])}
group<-extract('group');age<-as.numeric(extract('age'));gender<-extract('gender');ids<-sub('^(s[0-9]+).*','\\1',title);stopifnot(length(ids)==80,!anyDuplicated(ids))
counts<-read.delim(gzfile('work/geo_muscle/GSE144304_counts.txt.gz'),header=TRUE,check.names=FALSE,row.names=1)
sourceid<-sub('^(s[0-9]+).*','\\1',names(counts));idx<-match(sourceid,ids);stopifnot(!anyNA(idx),length(unique(idx))==80)
group<-group[idx];age<-age[idx];gender<-gender[idx];keep<-group %in% c('Fit','Frail') & gender=='Male'
stopifnot(sum(keep)==28,sum(group[keep]=='Frail')==13,sum(group[keep]=='Fit')==15,all(is.finite(age)),all(gender %in% c('Male','Female')),all(as.matrix(counts)>=0))
raw<-as.matrix(counts[,keep]);rownames(raw)<-sub('\\.[0-9]+$','',rownames(raw));stopifnot(!anyDuplicated(rownames(raw)));storage.mode(raw)<-'numeric'
design<-model.matrix(~I(group[keep]=='Frail')+age[keep]);colnames(design)<-c('intercept','frail','age')
y<-DGEList(counts=raw);ok<-filterByExpr(y,design=design);y<-y[ok,,keep.lib.sizes=FALSE];y<-calcNormFactors(y,method='TMM');y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);res<-glmQLFTest(fit,coef='frail');tt<-res$table
ids<-c(FOXP2='ENSG00000128573',SMAD5='ENSG00000113658',RBM6='ENSG00000004534',MEIS1='ENSG00000143995');stopifnot(all(ids%in%rownames(raw)))
df<-data.frame(gene=names(ids),ensembl=unname(ids),retained=ids %in% rownames(tt),logFC=NA_real_,logCPM=NA_real_,F=NA_real_,p=NA_real_,row.names=NULL)
for(i in seq_along(ids)) if(df$retained[i])df[i,c('logFC','logCPM','F','p')]<-as.numeric(tt[ids[i],c('logFC','logCPM','F','PValue')])
df$bh_p_4<-p.adjust(df$p,method='BH',n=4)
write.csv(df,'outputs/stage6_gse144304_edger_male_sensitivity.csv',row.names=FALSE,na='')
cat('n:',sum(keep),'input rows:',nrow(raw),'retained:',sum(ok),'\n');print(df,row.names=FALSE);cat('edgeR:',as.character(packageVersion('edgeR')),'\n')

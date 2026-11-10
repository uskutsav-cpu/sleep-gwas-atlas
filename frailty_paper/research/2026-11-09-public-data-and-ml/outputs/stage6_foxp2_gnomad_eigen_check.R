for(name in c('left','right')){
 d<-read.csv(paste0('outputs/stage6_foxp2_gnomad_fin_',name,'_variants.csv'));p<-nrow(d);R<-matrix(readBin(paste0('work/stage6_foxp2_gnomad_fin_',name,'_ld.bin'),'double',n=p*p,size=8,endian='little'),p,p);e<-eigen(R,symmetric=TRUE,only.values=TRUE)$values;cat(name,'n',p,'min_eigen',min(e),'negative_count',sum(e < -1e-8),'max',max(e),'\n');flush.console()
}

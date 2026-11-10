d<-read.csv('outputs/stage6_foxp2_gnomad_fin_left_variants.csv');n<-nrow(d);R<-matrix(readBin('work/stage6_foxp2_gnomad_fin_left_ld.bin','double',n=n*n,size=8,endian='little'),n,n)
for (radius in c(500000,400000,250000,100000)) {
 sel<-which(d$pos37>=114043159-radius & d$pos37<=114043159+radius);e<-eigen(R[sel,sel],symmetric=TRUE,only.values=TRUE)$values;cat('radius',radius,'variants',length(sel),'minimum_eigen',min(e),'negative_eigen',sum(e < -1e-8),'\n');flush.console()
}

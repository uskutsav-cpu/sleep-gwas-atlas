source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
d <- read.delim('outputs/shared-genetics-atlas/frailty_coloc/matched.tsv')
o <- list();i<-0
for(ch in c(3,7)) for(radius in c(100000,250000,500000)) {
 center <- if(ch==3) 117642005 else 114071035
 x <- d[d$chr==ch & abs(d$bp37-center)<=radius,]
 if(nrow(x)<500) {i<-i+1;o[[i]]<-data.frame(chr=ch,radius=radius,n=nrow(x),p12=NA,H3=NA,H4=NA);next}
 l1<-approx.bf.estimates(x$fi_beta/x$fi_se,x$fi_se^2,'quant')$lABF
 l2<-approx.bf.estimates(x$r13_beta/x$r13_se,x$r13_se^2,'cc')$lABF
 for(p12 in c(1e-5,1e-6)){z<-combine.abf(l1,l2,1e-4,1e-4,p12,quiet=TRUE);i<-i+1;o[[i]]<-data.frame(chr=ch,radius=radius,n=nrow(x),p12=p12,H3=unname(z[4]),H4=unname(z[5]))}
}
res<-do.call(rbind,o);write.table(res,'outputs/shared-genetics-atlas/frailty_coloc/window_sensitivity.tsv',sep='\t',quote=FALSE,row.names=FALSE);print(res)

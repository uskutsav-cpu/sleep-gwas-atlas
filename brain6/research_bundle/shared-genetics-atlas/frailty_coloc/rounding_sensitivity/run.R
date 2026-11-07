# Post-hoc stress test, not a probabilistic model of true rounding errors.
source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p <- 'outputs/shared-genetics-atlas/frailty_coloc/rounding_sensitivity'
a <- read.delim('outputs/shared-genetics-atlas/frailty_coloc/matched.tsv')
b <- read.delim('outputs/shared-genetics-atlas/frailty_coloc_chr5/matched.tsv')
d <- rbind(a,b); set.seed(20261027); out <- list(); k <- 0
for (chr in c(3,5,7)) {
 x <- d[d$chr==chr,]; stopifnot(nrow(x)>500, all(x$fi_se > 0.00005))
 r13 <- approx.bf.estimates(x$r13_beta/x$r13_se,x$r13_se^2,'cc')$lABF
 for (rep in 0:500) {
  # Uniform jitter is an illustrative perturbation, not an inferred measurement distribution.
  beta <- x$fi_beta; se <- x$fi_se
  if(rep>0) {beta <- beta+runif(nrow(x),-0.00005,0.00005); se <- se+runif(nrow(x),-0.00005,0.00005)}
  fi <- approx.bf.estimates(beta/se,se^2,'quant')$lABF
  lead <- x$rsid[which.max(fi+r13)]
  for (p12 in c(1e-5,1e-6)) {
   z <- combine.abf(fi,r13,1e-4,1e-4,p12,quiet=TRUE)
   k<-k+1;out[[k]]<-data.frame(chr=chr,n=nrow(x),iteration=rep,p12=p12,H3=unname(z[4]),H4=unname(z[5]),conditional_lead=lead)
  }
 }
}
res <- do.call(rbind,out)
write.table(res,file.path(p,'replicates.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
for(chr in c(3,5,7)) for(prior in c(1e-5,1e-6)) {
 x<-res[res$chr==chr & res$p12==prior,]; b<-x[x$iteration==0,]; y<-x[x$iteration>0,]
 cat('chr',chr,'p12',prior,'n',b$n,'baseline H4',b$H4,'jitter H4 q0/2.5/50/97.5/100',quantile(y$H4,c(0,.025,.5,.975,1)), 'lead baseline',b$conditional_lead,'unchanged',sum(y$conditional_lead==b$conditional_lead),'/500\n')
}

source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
root <- 'outputs/shared-genetics-atlas'; out <- file.path(root,'frailty_coloc/adhd_same_support')
f <- read.delim(file.path(root,'frailty_coloc/matched.tsv'), stringsAsFactors=FALSE)
f <- f[f$chr==3,]; p <- read.delim(file.path(root,'variant_context/rs62264764_dense_qtl/multisignal_pilot/variants.tsv'), stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(f$rsid), !anyDuplicated(p$rsid))
x <- merge(f,p,by='rsid',suffixes=c('_fi','_pilot'))
x <- x[x$bp37==x$pos37,]
direct <- x$effect_allele==x$alt & x$other_allele==x$ref
reverse <- x$effect_allele==x$ref & x$other_allele==x$alt
x <- x[direct|reverse,]; direct <- direct[direct|reverse]
x$adhd_aligned <- x$adhd_beta * ifelse(direct,1,-1)
x$brain_insomnia_aligned <- x$insomnia_beta * ifelse(direct,1,-1)
cols <- c('fi_beta','fi_se','r13_beta','r13_se','adhd_aligned','adhd_se','brain_insomnia_aligned','insomnia_se')
x <- x[apply(x[,cols],1,function(y) all(is.finite(y))) & x$fi_se>0 & x$r13_se>0 & x$adhd_se>0 & x$insomnia_se>0,]
stopifnot(nrow(x)>=500, !anyDuplicated(x$rsid), all(nchar(x$ref)==1),all(nchar(x$alt)==1))
fi <- approx.bf.estimates(x$fi_beta/x$fi_se,x$fi_se^2,'quant')$lABF
traits <- list(ADHD=c('adhd_aligned','adhd_se'), FinnGen_insomnia=c('r13_beta','r13_se'), Brain6_insomnia=c('brain_insomnia_aligned','insomnia_se'))
res <- list(); k <- 0
for (name in names(traits)) {
 b <- x[[traits[[name]][1]]]; se <- x[[traits[[name]][2]]]; second <- approx.bf.estimates(b/se,se^2,'cc')$lABF
 for(shared in c(1e-5,1e-6)) { z <- combine.abf(fi,second,1e-4,1e-4,shared,quiet=TRUE);stopifnot(abs(sum(z)-1)<1e-12);k<-k+1;res[[k]]<-data.frame(trait=name,n=nrow(x),p12=shared,H0=z[1],H1=z[2],H2=z[3],H3=z[4],H4=z[5],fi_min_wald_p=min(2*pnorm(-abs(x$fi_beta/x$fi_se))),trait_min_wald_p=min(2*pnorm(-abs(b/se)))) }
}
write.table(do.call(rbind,res),file.path(out,'posteriors.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('FI chr3 input',nrow(f),'pilot input',nrow(p),'intersection',nrow(x),'direct',sum(direct),'\n')
print(do.call(rbind,res))

source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p <- 'outputs/shared-genetics-atlas/frailty_coloc_chr5'
d <- read.delim(file.path(p,'matched.tsv'))
fits <- list(); top <- list(); i <- 0
for (ch in c(5)) {
 x <- d[d$chr==ch,]; if(nrow(x)<500) next
 for (prior in c('default','narrow','wide')) {
  if(prior=='default') { l1 <- approx.bf.estimates(x$fi_beta/x$fi_se,x$fi_se^2,'quant')$lABF; l2 <- approx.bf.estimates(x$r13_beta/x$r13_se,x$r13_se^2,'cc')$lABF }
  else {sd <- if(prior=='narrow') .1 else .3; calc <- function(beta,se) {r <- sd^2/(sd^2+se^2); .5*(log1p(-r)+r*(beta/se)^2)}; l1 <- calc(x$fi_beta,x$fi_se); l2 <- calc(x$r13_beta,x$r13_se)}
  for (p12 in c(1e-5,1e-6)) {
   z <- combine.abf(l1,l2,1e-4,1e-4,p12,quiet=TRUE); stopifnot(abs(sum(z)-1)<1e-12)
   i <- i+1; fits[[i]] <- data.frame(chr=ch,n=nrow(x),effect_prior=prior,p12=p12,H0=z[1],H1=z[2],H2=z[3],H3=z[4],H4=z[5],row.names=NULL)
  }
  if(prior=='default') { x$conditional_h4 <- exp(l1+l2-logsum(l1+l2)); x <- x[order(-x$conditional_h4),]; x$credible95 <- cumsum(x$conditional_h4)-x$conditional_h4 < .95; top[[as.character(ch)]] <- x }
 }
}
write.table(do.call(rbind,fits),file.path(p,'posteriors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,top),file.path(p,'conditional_variant_posteriors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(do.call(rbind,fits));for(ch in names(top)) { cat('\nchr',ch,'top\n');print(head(top[[ch]][,c('rsid','bp37','conditional_h4','fi_p','r13_p')],5));cat('95set size',sum(top[[ch]]$credible95),'\n') }

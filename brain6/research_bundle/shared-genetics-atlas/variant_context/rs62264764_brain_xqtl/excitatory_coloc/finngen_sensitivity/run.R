source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_brain_xqtl/excitatory_coloc/finngen_sensitivity'
d <- read.delim(file.path(p,'matched.tsv'))
stopifnot(nrow(d)>=500, !anyDuplicated(d$variant), all(d$qtl_se>0), all(d$gwas_se>0), all(d$fg_se>0))
x <- approx.bf.estimates(d$qtl_beta/d$qtl_se,d$qtl_se^2,'quant')$lABF
res <- list()
for(trait in c('brain6','finngen')) {
 b <- if(trait=='brain6') d$gwas_beta_alt else d$fg_beta_alt
 se <- if(trait=='brain6') d$gwas_se else d$fg_se
 y <- approx.bf.estimates(b/se,se^2,'cc')$lABF
 for(p12 in c(1e-5,1e-6)) {
  v <- combine.abf(x,y,1e-4,1e-4,p12,quiet=TRUE)
  stopifnot(abs(sum(v)-1)<1e-10)
  res[[length(res)+1]] <- data.frame(trait=trait,n=nrow(d),p12=p12,H0=v[1],H1=v[2],H2=v[3],H3=v[4],H4=v[5],row.names=NULL)
 }
}
results <- do.call(rbind,res)
write.table(results,file.path(p,'posteriors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(results)

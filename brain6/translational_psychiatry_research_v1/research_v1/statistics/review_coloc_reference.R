# Execute unchanged retained coloc 5.2.3 reference functions on the existing sensitivity grid.
# This is cross-language numerical validation, not a full native coloc package rerun.
options(digits=22)
args <- commandArgs(trailingOnly=TRUE)
repo <- normalizePath(args[1],mustWork=TRUE)
base <- file.path(repo,'brain6/translational_psychiatry_research_v1/research_v1')
source(file.path(base,'molecular/coloc_5.2.3_source.R'))
conditions <- read.delim(file.path(base,'molecular/sensitivity_results.tsv'),check.names=FALSE)
rows <- list(); idx <- 0
for (path in unique(conditions$input_path)) {
  d <- read.delim(gzfile(file.path(repo,path)),check.names=FALSE)
  sd <- sdY.est(d$qtl_se^2,d$qtl_maf,d$qtl_n)
  g <- approx.bf.estimates(z=d$gwas_beta/d$gwas_se,V=d$gwas_se^2,type='cc')$lABF
  for (scale in c(.5,1,2)) {
    q <- approx.bf.estimates(z=d$qtl_beta/d$qtl_se,V=d$qtl_se^2,type='quant',sdY=sd*scale)$lABF
    for (p12 in c(1e-6,5e-6,1e-5,5e-5,1e-4)) {
      pp <- combine.abf(g,q,1e-4,1e-4,p12,quiet=TRUE)
      old <- conditions[conditions$input_path==path & conditions$sdY_scale==scale & conditions$p12==p12,]
      stopifnot(nrow(old)==1)
      err <- max(abs(pp-as.numeric(old[paste0('pp_h',0:4)])))
      stopifnot(err<1e-10)
      idx <- idx+1
      rows[[idx]] <- data.frame(input_path=path,sdY_scale=scale,p12=p12,reference_R_sdY=sd,reference_R_H4=pp[5],max_posterior_absolute_error=err,classification='REFERENCE_R_FUNCTION_NUMERICAL_EXECUTION_NOT_FULL_PACKAGE',stringsAsFactors=FALSE)
    }
  }
}
out <- do.call(rbind,rows)
stopifnot(nrow(out)==60)
write.table(out,file.path(base,'statistics/reference_R_posterior_checks.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(base,'statistics/reference_R_environment.txt'))
cat('Reference R functions:',nrow(out),'conditions;maximum error',max(out$max_posterior_absolute_error),'\n')

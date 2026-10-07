# Read saved fitted objects and actual signed LD; no fitting and no source-row export.
options(digits=22)
repo<-normalizePath(commandArgs(trailingOnly=TRUE)[1],mustWork=TRUE)
stat<-file.path(repo,'brain6/translational_psychiatry_research_v1/research_v1/statistics')
ldout<-file.path(repo,'brain6/translational_psychiatry_research_v1/research_v1/ld')
local<-file.path(repo,'work/ld_genotypes_research_v1/C_finemap_inputs')
summaries<-read.delim(file.path(ldout,'native_finemap_fit_summary.tsv'),check.names=FALSE)
cred<-read.delim(file.path(ldout,'native_finemap_credible_sets.tsv'),check.names=FALSE)
pips<-read.delim(file.path(ldout,'native_finemap_variant_PIPs.tsv'),check.names=FALSE)
history<-read.delim(file.path(ldout,'native_finemap_ELBO_history.tsv'),check.names=FALSE)
p<-unique(summaries$n_snps);stopifnot(length(p)==1,p==2185,nrow(summaries)==8)
R<-matrix(readBin(file.path(local,'shared_signed_LD.float64.bin'),what=numeric(),n=p*p,size=8,endian='little'),nrow=p)
stopifnot(all(is.finite(R)),max(abs(R-t(R)))<1e-12,max(abs(diag(R)-1))<1e-12,max(abs(R))<=1+1e-12)
checks<-list();idx<-0
for(i in seq_len(nrow(summaries))){
  row<-summaries[i,];fit<-readRDS(file.path(local,paste0(row$trait,'_',row$model,'_native_fit.rds')))
  stopifnot(isTRUE(fit$converged),all(is.finite(fit$alpha)),all(fit$alpha>=0&fit$alpha<=1),max(abs(rowSums(fit$alpha)-1))<1e-12)
  independent_pip<-1-exp(colSums(log1p(-fit$alpha)))
  pi<-pips[pips$model==row$model&pips$trait==row$trait,]
  stopifnot(nrow(pi)==p,identical(as.character(pi$variant_key),colnames(fit$alpha)))
  perr<-max(abs(independent_pip-pi$PIP))
  hh<-history[history$model==row$model&history$trait==row$trait,]
  eerr<-max(abs(hh$ELBO-fit$elbo));minchange<-if(length(fit$elbo)>1)min(diff(fit$elbo)) else 0
  stopifnot(perr<1e-12,eerr<1e-7,minchange>=-1e-6,length(fit$sets$cs)==row$n_credible_sets)
  for(j in seq_along(fit$sets$cs)){
    cs<-fit$sets$cs[[j]];effect<-fit$sets$cs_index[j];csname<-names(fit$sets$cs)[j]
    coverage<-sum(fit$alpha[effect,cs])
    sorted<-sort(fit$alpha[effect,],decreasing=TRUE)
    minsize<-which(cumsum(sorted)>=.95)[1]
    sub<-abs(R[cs,cs,drop=FALSE]);values<-sub[upper.tri(sub)]
    if(length(values)==0)values<-1
    purity<-c(min(values),mean(values),median(values))
    tt<-cred[cred$model==row$model&cred$trait==row$trait&cred$credible_set==csname,]
    stopifnot(nrow(tt)==length(cs),setequal(tt$variant_key,colnames(fit$alpha)[cs]))
    tidx<-match(tt$variant_key,colnames(fit$alpha))
    alpha_error<-max(abs(tt$within_signal_posterior-fit$alpha[effect,tidx]))
    purity_error<-max(abs(purity-c(tt$purity_min_abs_corr[1],tt$purity_mean_abs_corr[1],tt$purity_median_abs_corr[1])))
    stopifnot(coverage>=.95-1e-12,length(cs)==minsize,purity[1]>=.5,alpha_error<1e-12,purity_error<1e-12)
    idx<-idx+1
    checks[[idx]]<-data.frame(model=row$model,trait=row$trait,credible_set=csname,effect_index=effect,CS_size=length(cs),independent_conditional_coverage=coverage,native_conditional_coverage=fit$sets$coverage[j],independent_purity_min=purity[1],independent_purity_mean=purity[2],independent_purity_median=purity[3],max_PIP_error=perr,max_CS_alpha_error=alpha_error,max_purity_error=purity_error,max_ELBO_roundtrip_error=eerr,min_ELBO_change=minchange,alpha_row_sum_max_error=max(abs(rowSums(fit$alpha)-1)),classification='POST_RESULT_READ_ONLY_SAVED_FIT_AND_LD_REVIEW',independent_cohort_replication=FALSE)
  }
}
stopifnot(length(checks)==8)
write.table(do.call(rbind,checks),file.path(stat,'native_fine_mapping_RDS_checks.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(stat,'native_fine_mapping_review_R_environment.txt'))
cat('Independently reviewed8savedfits,8qualifyingcredible sets,actualLDpurity,PIPandELBO;no fits rerun\n')

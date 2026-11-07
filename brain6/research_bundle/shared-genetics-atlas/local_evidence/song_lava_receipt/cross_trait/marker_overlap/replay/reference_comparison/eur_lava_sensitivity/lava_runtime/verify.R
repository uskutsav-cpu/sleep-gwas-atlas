suppressPackageStartupMessages(library(jsonlite))
p<-'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay/reference_comparison/eur_lava_sensitivity'
z<-fromJSON(file.path(p,'lava_runtime/lava_result.json'));e<-fromJSON(file.path(p,'eur_result.json'))
stopifnot(z$n_snps==e$matched_markers,z$K==e$K_99,identical(z$failed,c(FALSE,FALSE)))
stopifnot(abs(z$omega[1,1]-e$fi_omega)<1e-12)
stat<-z$delta_squared_sums / diag(z$sigma)
pval<-c(pf(stat[1]/z$K,z$K,z$N[1]-z$K-1,lower.tail=FALSE),pchisq(stat[2],z$K,lower.tail=FALSE))
stopifnot(all(abs((pval-z$univ$p)/z$univ$p)<1e-5))
cat('Verified two univariate p-values and prior FI projection. This does not validate reference calibration.\n')

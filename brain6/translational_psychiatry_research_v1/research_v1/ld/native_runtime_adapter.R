# Pinned susieR0.14.2 interface; never silently treat absent diagnostics as pass.
extract_kriging_table <- function(result,p) {
 tab <- result$conditional_dist
 required <- c('z','condmean','condvar','z_std_diff','logLR')
 if(!is.data.frame(tab)||nrow(tab)!=p||!all(required%in%names(tab)))stop('Invalid pinned kriging_rss conditional_dist interface')
 if(!all(is.finite(as.matrix(tab))))stop('Nonfinite native kriging diagnostics')
 tab
}
rss_large_n_consistency <- function(z,R) {
 # In0.14.2 explicit n=Inf yields Inf/Inf; omitting n implements its documented large-N limit.
 susieR::estimate_s_rss(z=z,R=R,method='null-mle')
}

p <- 'outputs/shared-genetics-atlas/global_audit/brain6_frailty_covariance/three_frailty_common_support'
x <- read.delim(file.path(p,'delete_rg.tsv')); y <- read.delim(file.path(p,'contrasts.tsv')); stopifnot(nrow(x)==200,nrow(y)==3)
for (j in seq_len(nrow(y))) {
 tr <- strsplit(y$contrast[j], ' minus ')[[1]]; d <- x[[tr[1]]] - x[[tr[2]]]
 se <- sqrt((length(d)-1)/length(d)*sum((d-mean(d))^2))
 stopifnot(abs(se-y$se_paired[j])<1e-12,abs(min(1,3*2*pnorm(-abs(y$difference[j]/se)))-y$bonferroni_p_3[j])<1e-12)
}
cat('R independently reproduces all three paired SEs and Bonferroni p-values.\n')

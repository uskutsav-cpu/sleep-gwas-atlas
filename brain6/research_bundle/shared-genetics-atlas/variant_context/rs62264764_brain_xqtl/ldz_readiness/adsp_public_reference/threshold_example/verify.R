p <- 'outputs/shared-genetics-atlas/variant_context/rs62264764_brain_xqtl/ldz_readiness/adsp_public_reference/threshold_example'
R <- matrix(c(1,.81,.60,.81,1,.19,.60,.19,1),nrow=3,byrow=TRUE)
Z <- R; Z[2,3] <- 0; Z[3,2] <- 0
r <- data.frame(matrix=c('unthresholded','missing_pair_zero_filled'),min_eigenvalue=c(min(eigen(R,symmetric=TRUE)$values),min(eigen(Z,symmetric=TRUE)$values)),determinant=c(det(R),det(Z)))
stopifnot(r$min_eigenvalue[1]>0,r$min_eigenvalue[2]<0)
write.table(r,file.path(p,'results.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(r)

p<-'outputs/shared-genetics-atlas/variant_context/rs62264764_dense_qtl/multisignal_pilot'
R<-as.matrix(read.delim(gzfile(file.path(p,'ld_matrix.tsv.gz')),header=FALSE));e<-eigen(R,symmetric=TRUE,only.values=TRUE)$values
out<-data.frame(n=length(e),min_eigenvalue=min(e),max_eigenvalue=max(e),rank_at_1e_minus_8=sum(e>1e-8),eigenvalues_below_minus_1e_minus_8=sum(e< -1e-8))
write.table(out,file.path(p,'ld_eigen_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE);print(out)

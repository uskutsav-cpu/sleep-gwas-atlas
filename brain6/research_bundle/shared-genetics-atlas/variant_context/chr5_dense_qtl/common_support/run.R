source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p <- 'outputs/shared-genetics-atlas/variant_context/chr5_dense_qtl'; out <- file.path(p,'common_support')
a<-read.delim(file.path(p,'Testis_ENSG00000251574_insomnia_matched.tsv'),stringsAsFactors=FALSE)
b<-read.delim(file.path(p,'Testis_ENSG00000251574_adhd_matched.tsv'),stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(a$variant),!anyDuplicated(b$variant))
x<-merge(a,b,by='variant',suffixes=c('_i','_a'))
stopifnot(nrow(x)>=500,all(x$ref_i==x$ref_a),all(x$alt_i==x$alt_a),all(x$qtl_beta_i==x$qtl_beta_a),all(x$qtl_se_i==x$qtl_se_a))
score<-function(beta,se,type) approx.bf.estimates(beta/se,se^2,type)$lABF
q<-score(x$qtl_beta_i,x$qtl_se_i,'quant');i<-score(x$gwas_beta_alt_i,x$gwas_se_i,'cc');a<-score(x$gwas_beta_alt_a,x$gwas_se_a,'cc')
pairs<-list(qtl_insomnia=list(q,i),qtl_adhd=list(q,a),insomnia_adhd=list(i,a));r<-list();credible<-list();k<-0
for(name in names(pairs)){
 s<-pairs[[name]][[1]];t<-pairs[[name]][[2]]
 for(p12 in c(1e-5,1e-6)){z<-combine.abf(s,t,1e-4,1e-4,p12,quiet=TRUE);stopifnot(abs(sum(z)-1)<1e-12);k<-k+1;r[[k]]<-data.frame(pair=name,n=nrow(x),p12=p12,H0=z[1],H1=z[2],H2=z[3],H3=z[4],H4=z[5],row.names=NULL)}
 posterior<-exp(s+t-logsum(s+t));order_idx<-order(-posterior,x$variant);inc<-cumsum(posterior[order_idx])-posterior[order_idx]<.95
 y<-data.frame(pair=name,variant=x$variant[order_idx],rsid=x$qtl_rsid_i[order_idx],conditional_h4=posterior[order_idx],credible95=inc)
 credible[[name]]<-y
}
write.table(do.call(rbind,r),file.path(out,'pair_posteriors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,credible),file.path(out,'conditional_snps.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
sets<-lapply(credible,function(z)z$variant[z$credible95]);cat('Common variants:',nrow(x),'\n');print(do.call(rbind,r));for(name in names(sets)) {cat(name,'95-set',length(sets[[name]]),'leader',credible[[name]]$rsid[1],'prob',credible[[name]]$conditional_h4[1],'\n')};cat('all-three intersection',length(Reduce(intersect,sets)),'\n')

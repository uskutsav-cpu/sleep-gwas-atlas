source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p <- 'outputs/shared-genetics-atlas/component_exclusion/schoeler_physical_function/chr5_coloc'
i <- read.delim(file.path(p,'insomnia_matched.tsv'));a <- read.delim(file.path(p,'adhd_matched.tsv'))
stopifnot(!anyDuplicated(i$pos38),!anyDuplicated(a$pos38))
i<-i[order(i$pos38),];a<-a[order(a$pos38),];m<-intersect(i$pos38,a$pos38)
i<-i[match(m,i$pos38),];a<-a[match(m,a$pos38),]
stopifnot(nrow(i)==length(m),nrow(a)==length(m),length(m)>=500,identical(i$effect_allele,a$effect_allele),identical(i$other_allele,a$other_allele),isTRUE(all.equal(i$physical_beta,a$physical_beta)),isTRUE(all.equal(i$physical_se,a$physical_se)))
result<-list();j<-0
for(trait in c('insomnia','adhd')) {
 d<-if(trait=='insomnia')i else a
 x<-approx.bf.estimates(d$physical_beta/d$physical_se,d$physical_se^2,'quant')$lABF
 y<-approx.bf.estimates(d$trait_beta/d$trait_se,d$trait_se^2,'cc')$lABF
 for(q in c(1e-5,1e-6)) {
  z<-combine.abf(x,y,1e-4,1e-4,q,quiet=TRUE);stopifnot(abs(sum(z)-1)<1e-12)
  j<-j+1;result[[j]]<-data.frame(trait=trait,n=nrow(d),p12=q,H0=z[1],H1=z[2],H2=z[3],H3=z[4],H4=z[5],row.names=NULL)
 }
}
res<-do.call(rbind,result);write.table(res,file.path(p,'common_support/posteriors.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(m,file.path(p,'common_support/positions.tsv'),sep='\t',row.names=FALSE,quote=FALSE,col.names='pos38')
print(res)

source('outputs/shared-genetics-atlas/frailty_coloc/coloc_5.2.3_claudia.R')
p<-'outputs/shared-genetics-atlas/variant_context/rs62264764_brain_xqtl/dense_coloc'
base<-'/Volumes/Brain6Work/xqtl-chr3/matched'
res<-list();tops<-list();i<-0L
for(ctx in c('oli'))for(trait in c('insomnia','adhd')){
 d<-read.delim(file.path(base,paste0(ctx,'_',trait,'.tsv')));stopifnot(nrow(d)>=500,!anyDuplicated(d$variant),all(d$qtl_se>0),all(d$gwas_se>0))
 l1<-approx.bf.estimates(d$qtl_beta/d$qtl_se,d$qtl_se^2,'quant')$lABF
 l2<-approx.bf.estimates(d$gwas_beta_alt/d$gwas_se,d$gwas_se^2,'cc')$lABF
 for(p12 in c(1e-5,1e-6)){
  v<-combine.abf(l1,l2,1e-4,1e-4,p12,quiet=TRUE);stopifnot(abs(sum(v)-1)<1e-12);i<-i+1L
  res[[i]]<-data.frame(cell=ctx,trait=trait,n=nrow(d),p12=p12,H0=v[1],H1=v[2],H2=v[3],H3=v[4],H4=v[5],row.names=NULL)
 }
 cs<-exp(l1+l2-logsum(l1+l2));ord<-order(-cs);sel<-ord[seq_len(min(10,length(ord)))];tops[[length(tops)+1L]]<-data.frame(cell=ctx,trait=trait,variant=d$variant[sel],rsid=d$gwas_rsid[sel],conditional_h4=cs[sel])
}
write.table(do.call(rbind,res),file.path(p,'posteriors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,tops),file.path(p,'top_conditional.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(do.call(rbind,res));writeLines(capture.output(sessionInfo()),file.path(p,'session_info.txt'))

p <- 'outputs/shared-genetics-atlas/conditional_enrichment'
root <- '/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support'
s <- read.delim(gzfile(file.path(root,'S.sumstats.gz')),stringsAsFactors=FALSE)
set.seed(20261027); n_sim<-1000000L; out<-list(); i<-0L
pars<-list(FI=c(1.271851342700776,1.4045250827610078,.08081972497537657),FFS=c(1.2714481467007943,1.5006275551853958,.052417344112559244))
for(pair in c('FI','FFS')) {
 t<-read.delim(gzfile(file.path(root,paste0(pair,'.sumstats.gz'))),stringsAsFactors=FALSE)
 stopifnot(nrow(s)==1138262,nrow(t)==nrow(s),!anyDuplicated(s$SNP),identical(s$SNP,t$SNP),identical(s$A1,t$A1),identical(s$A2,t$A2),all(is.finite(s$Z)),all(is.finite(t$Z)))
 rho<-pars[[pair]][3]; z1<-rnorm(n_sim); z2<-rho*z1+sqrt(1-rho^2)*rnorm(n_sim)
 models<-list(observed=list(s$Z,t$Z),unit_null=list(z1,z2),inflated_zero_mean=list(z1*sqrt(pars[[pair]][1]),z2*sqrt(pars[[pair]][2])))
 for(model in names(models)) for(direction in c('insomnia_given_frailty','frailty_given_insomnia')) {
  z<-models[[model]]; ps<-2*pnorm(-abs(z[[1]])); pt<-2*pnorm(-abs(z[[2]]));target<-if(direction=='insomnia_given_frailty')ps else pt; condition<-if(direction=='insomnia_given_frailty')pt else ps
  for(cut in c(1,.1,.01,.001)) for(target_cut in c(.05,.01,.001)) {
   in_subset<-condition<=cut; denom<-sum(in_subset);num<-sum(in_subset & target<=target_cut);i<-i+1L
   out[[i]]<-data.frame(pair=pair,model=model,direction=direction,condition_p_max=cut,target_p_max=target_cut,selected=denom,target_count=num,fraction=if(denom)num/denom else NA_real_)
  }
 }
}
d<-do.call(rbind,out);write.table(d,file.path(p,'conditional_fractions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('rows',nrow(d),'observed support',nrow(s),'simulations per pair per model',n_sim,'\n')

#!/usr/bin/env Rscript
out <- 'outputs/shared-genetics-atlas/conditional_enrichment/conjfdr_readiness/r_implementation/synthetic_ld_overlap/maxp_boundary_null'
n <- 44000L; nb <- 4400L; block <- rep(seq_len(nb),each=10L); rows <- list()
for (seed in 82001:83000) {
 set.seed(seed)
 e1 <- sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n)
 e2 <- sqrt(.5)*rnorm(nb)[block]+sqrt(.5)*rnorm(n)
 p1 <- 2*pnorm(-abs(e1+12))
 for (rho in c(0,.3)) {
  p2 <- 2*pnorm(-abs(rho*e1+sqrt(1-rho^2)*e2)); p <- pmax(p1,p2)
  stopifnot(all(p>=p1),all(p>=p2),all(p>=0 & p<=1))
  for(method in c('bonferroni','BY','BH')) {
   adj<-p.adjust(p,method=method);stopifnot(all(adj>=0 & adj<=1))
   if(seed==82001L){s<-sort(p);k<-seq_len(n);expected<-switch(method,bonferroni=sum(p<=.05/n),BH=max(c(0L,k[s<=.05*k/n])),BY=max(c(0L,k[s<=.05*k/(n*sum(1/k))])));stopifnot(sum(adj<=.05)==expected)}
   rows[[length(rows)+1L]]<-data.frame(seed=seed,rho=rho,method=method,rejections=sum(adj<=.05),minimum=min(adj))
  }
 }
 if(seed%%100==0){cat('completed',seed,'\n');flush.console()}
}
d<-do.call(rbind,rows);stopifnot(nrow(d)==6000L)
write.table(d,file.path(out,'runs.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
summary<-do.call(rbind,lapply(split(d,list(d$rho,d$method),drop=TRUE),function(x){events<-sum(x$rejections>0);ci<-binom.test(events,nrow(x))$conf.int;data.frame(rho=x$rho[1],method=x$method[1],runs=nrow(x),events=events,rate=events/nrow(x),low=ci[1],high=ci[2],max_rejections=max(x$rejections))}))
write.table(summary,file.path(out,'summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS: direct first-seed counts agree; 6000 records saved.\n');print(sessionInfo())

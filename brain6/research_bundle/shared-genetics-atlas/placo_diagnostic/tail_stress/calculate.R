args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1);out<-args[1]
source('/Volumes/Extreme SSD/brain6-work/PLACO-v0.2.0/PLACO_v0.2.0.R')
params<-list(FI=c(1.271851342700776,1.4045250827610078,0.08081972497537657),FFS=c(1.2714481467007943,1.5006275551853958,0.052417344112559244))
alpha<-c(nominal=.05,screen=.05/(2*1138262));cap<-sqrt(80)
placop<-function(t,p) placo.plus(c(t,1),VarZ=p[1:2],CorZ=p[3],AbsTol=1e-14)$p.placo.plus
integral<-function(t,mu,rho,tol=1e-15,alternate=FALSE) {
 if(t>=cap^2)return(c(prob=0,err=0))
 a<-t/cap; sd<-sqrt(1-rho^2)
 f<-function(x) {
  ay<-t/abs(x); cy<-rho*(x-mu)
  # P(ay <= Y <= cap) + P(-cap <= Y <= -ay), conditional on X=x.
  pp<-pnorm((ay-cy)/sd,lower.tail=FALSE)-pnorm((cap-cy)/sd,lower.tail=FALSE)
  pm<-pnorm((-ay-cy)/sd)-pnorm((-cap-cy)/sd)
  dnorm(x,mean=mu)*(pmax(0,pp)+pmax(0,pm))
 }
 if(alternate) {g<-function(u) f(u+mu); lo<-c(-cap-mu,a-mu);hi<-c(-a-mu,cap-mu)} else {g<-f;lo<-c(-cap,a);hi<-c(-a,cap)}
 z<-lapply(seq_along(lo),function(i)integrate(g,lo[i],hi[i],abs.tol=tol,rel.tol=1e-10,subdivisions=1000,stop.on.error=TRUE))
 c(prob=sum(vapply(z,`[[`,numeric(1),'value')),err=sum(vapply(z,`[[`,numeric(1),'abs.error')))
}
# Independent conditional integration over the null variable Y instead of X.
integral_y<-function(t,mu,rho) {
 if(t>=cap^2)return(0)
 a<-t/cap; sd<-sqrt(1-rho^2)
 f<-function(y) {
  ax<-t/abs(y); cx<-mu+rho*y
  pp<-pnorm((ax-cx)/sd,lower.tail=FALSE)-pnorm((cap-cx)/sd,lower.tail=FALSE)
  pm<-pnorm((-ax-cx)/sd)-pnorm((-cap-cx)/sd)
  dnorm(y)*(pmax(0,pp)+pmax(0,pm))
 }
 sum(vapply(list(c(-cap,-a),c(a,cap)),function(interval)integrate(f,interval[1],interval[2],abs.tol=1e-17,rel.tol=1e-10,subdivisions=1000)$value,numeric(1)))
}
res<-list();i<-0
for(pair in names(params))for(lev in names(alpha)) {
 par<-params[[pair]];target<-alpha[[lev]]
 t<-uniroot(function(x)placop(x,par)-target,interval=c(.05,79.9),tol=1e-10)$root
 check<-placop(t,par);stopifnot(abs(check-target)<max(2e-13,target*1e-5))
 for(mu in 0:6){z<-integral(t,mu,par[3]);alt<-integral(t,mu,par[3],tol=1e-17,alternate=TRUE);sym<-integral(t,-mu,par[3]);other<-integral_y(t,mu,par[3]);stopifnot(abs(z[1]-other)<max(1e-13,z[1]*1e-7));stopifnot(abs(z[1]-alt[1])<max(1e-13,z[1]*1e-7),abs(z[1]-sym[1])<max(1e-13,z[1]*1e-7));i<-i+1;res[[i]]<-data.frame(pair=pair,level=lev,alpha=target,critical_product=t,mean_nonnull=mu,probability=z[1],quadrature_abs_error=z[2],alternate_probability=alt[1],alternate_abs_error=alt[2],independent_condition_probability=other,ratio_to_alpha=z[1]/target)}
}
d<-do.call(rbind,res);write.table(d,file.path(out,'tail_probabilities.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
old<-read.delim('outputs/shared-genetics-atlas/placo_diagnostic/simulated_null.tsv');for(pr in names(params))for(sc in c('sleep_only','frailty_only')) {x<-old[old$pair==pr & old$scenario==sc,];q<-d$probability[d$pair==pr & d$level=='nominal' & d$mean_nonnull==2];obs<-mean(x$P<.05);stopifnot(nrow(x)==1000,abs(obs-q)<5*sqrt(q*(1-q)/1000)+.002);cat('simulation comparison',pr,sc,'observed',obs,'quadrature',q,'\n')}
print(d,row.names=FALSE)

args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1);out<-args[1]
source('/Volumes/Extreme SSD/brain6-work/PLACO-v0.2.0/PLACO_v0.2.0.R')
params<-list(FI=c(1.271851342700776,1.4045250827610078,0.08081972497537657),FFS=c(1.2714481467007943,1.5006275551853958,0.052417344112559244))
alphas<-c(nominal=.05,screen=.05/(2*1138262));cap<-sqrt(80)
tail_prob<-function(t,sx,sy,rho,capped,swap=FALSE){
 if(swap){tmp<-sx;sx<-sy;sy<-tmp}
 upper<-if(capped)cap else Inf;lower<-if(capped)t/cap else 0
 f<-function(x){
  lim<-t/abs(x);mean<-rho*sy/sx*x;sd<-sy*sqrt(1-rho*rho)
  positive<-pnorm((lim-mean)/sd,lower.tail=FALSE)-pnorm((upper-mean)/sd,lower.tail=FALSE)
  negative<-pnorm((-lim-mean)/sd)-pnorm((-upper-mean)/sd)
  dnorm(x,sd=sx)*pmax(0,positive+negative)
 }
 pieces<-list(c(-upper,-lower),c(lower,upper))
 sum(vapply(pieces,function(z)integrate(f,z[1],z[2],abs.tol=1e-18,rel.tol=1e-9,subdivisions=1000,stop.on.error=TRUE)$value,numeric(1)))
}
outrows<-list()
for(pair in names(params))for(level in names(alphas)){
 p<-params[[pair]];alpha<-alphas[[level]]
 f<-function(t)placo.plus(c(t,1),VarZ=p[1:2],CorZ=p[3],AbsTol=1e-15)$p.placo.plus
 threshold<-uniroot(function(t)f(t)-alpha,c(.05,79.9),tol=1e-11)$root
 specs<-list(H00=c(1,1),H01=c(1,sqrt(p[2])),H02=c(sqrt(p[1]),1));prob<-list()
 for(name in names(specs)){
  s<-specs[[name]];u<-tail_prob(threshold,s[1],s[2],p[3],FALSE);us<-tail_prob(threshold,s[1],s[2],p[3],FALSE,TRUE)
  v<-tail_prob(threshold,s[1],s[2],p[3],TRUE);vs<-tail_prob(threshold,s[1],s[2],p[3],TRUE,TRUE)
  stopifnot(abs(u-us)<max(1e-15,u*1e-7),abs(v-vs)<max(1e-15,v*1e-7),v<=u+1e-15)
  prob[[name]]<-u
  outrows[[length(outrows)+1]]<-data.frame(pair=pair,level=level,alpha=alpha,critical_product=threshold,subnull=name,uncapped_probability=u,capped_probability=v,uncapped_over_alpha=u/alpha,capped_over_alpha=v/alpha,alternate_uncapped=us,alternate_capped=vs)
 }
 reconstructed<-prob$H01+prob$H02-prob$H00
 stopifnot(abs(reconstructed-f(threshold))<max(1e-15,alpha*1e-6))
 cat(pair,level,'official=',f(threshold),'reconstructed=',reconstructed,'\n')
}
write.table(do.call(rbind,outrows),file.path(out,'model_probabilities.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'runtime.txt'))

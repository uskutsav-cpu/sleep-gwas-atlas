d<-read.csv(file('stdin'),na.strings=c('','NA'));r<-read.delim('outputs/shared-genetics-atlas/component_exclusion/nhanes_public/prepandemic_availability/physical_sleep_pilot/results.tsv');wc<-function(x,y,w){mx<-weighted.mean(x,w);my<-weighted.mean(y,w);sum(w*(x-mx)*(y-my))/sqrt(sum(w*(x-mx)^2)*sum(w*(y-my)^2))}
for(scenario in unique(r$scenario)){
 x<-d;if(scenario!='exclude_ambiguous')x$weightloss[x$ambiguous==1]<-if(scenario=='ambiguous_as_one')1 else 0
 ok<-x$RIDAGEYR>=60 & complete.cases(x[,c('sleep','weightloss','inactivity','energy')]);x<-x[ok,];a<-x$weightloss+x$inactivity;b<-a+x$energy;aa<-wc(x$sleep,a,x$WTMECPRP);bb<-wc(x$sleep,b,x$WTMECPRP);v<-r[r$scenario==scenario,];stopifnot(nrow(x)==v$n[1],max(abs(v$estimate-c(aa,bb,bb-aa)))<1e-10)
}
cat('Weighted correlations and sample counts independently reproduced for all scenarios.\n')

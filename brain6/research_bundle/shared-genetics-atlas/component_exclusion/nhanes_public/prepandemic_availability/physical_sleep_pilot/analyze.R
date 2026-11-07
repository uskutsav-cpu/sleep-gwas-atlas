suppressPackageStartupMessages(library(survey))
dat<-read.csv(file('stdin'),na.strings=c('','NA'));stopifnot(all(dat$WTMECPRP>0));out<-list()
for(scenario in c('exclude_ambiguous','ambiguous_as_zero','ambiguous_as_one')){
 d<-dat
 if(scenario!='exclude_ambiguous') d$weightloss[d$ambiguous==1]<-if(scenario=='ambiguous_as_one') 1 else 0
 d$a<-d$weightloss+d$inactivity;d$b<-d$a+d$energy;d$s<-d$sleep;d$s2<-d$s^2;d$a2<-d$a^2;d$b2<-d$b^2;d$sa<-d$s*d$a;d$sb<-d$s*d$b
 d$eligible<-d$RIDAGEYR>=60 & complete.cases(d[,c('s','a','b')]);des<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMECPRP,nest=TRUE,data=d);dd<-subset(des,eligible);m<-svymean(~s+s2+a+a2+b+b2+sa+sb,dd)
 ra<-quote((sa-s*a)/sqrt((s2-s^2)*(a2-a^2)));rb<-quote((sb-s*b)/sqrt((s2-s^2)*(b2-b^2)));dif<-substitute(B-A,list(B=rb,A=ra));q<-list(svycontrast(m,ra),svycontrast(m,rb),svycontrast(m,dif));df<-degf(dd);crit<-qt(.975,df)
 for(i in seq_along(q)){est<-unname(coef(q[[i]]));se<-unname(SE(q[[i]]));out[[length(out)+1]]<-data.frame(scenario=scenario,n=sum(d$eligible),design_df=df,quantity=c('two_item_r','three_item_r','three_minus_two')[i],estimate=est,se=se,lower=est-crit*se,upper=est+crit*se)}
}
write.table(do.call(rbind,out),'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/prepandemic_availability/physical_sleep_pilot/results.tsv',sep='\t',row.names=FALSE,quote=FALSE)

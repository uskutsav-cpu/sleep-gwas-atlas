#!/usr/bin/env Rscript
# Run from this directory with survey, survival installed. Independent implementation in R survey.
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='fail')
d <- read.csv('survey_full_design_frame.csv',na.strings=c('','NA'))
d$cycle_f <- factor(d$cycle,levels=c('2011-2012','2013-2014'));d$sex_f<-factor(d$RIAGENDR,levels=c(1,2),labels=c('male','female'));d$race_f<-factor(d$RIDRETH3);d$edu_f<-factor(d$education_clean);d$age_c<-d$RIDAGEYR-70;d$fi_diff_nosleep<-d$fi_all-d$fi_no_sleep;d$fi_diff_nofatigue<-d$fi_all-d$fi_no_fatigue;d$fi_10<-d$fi_neither*10;d$fi_all10<-d$fi_all*10;d$fi_nomoodcare10<-d$fi_no_mood_no_care*10
des<-svydesign(ids=~SDMVPSU,strata=~stratum,weights=~wt4,data=d,nest=TRUE);stopifnot(degf(des)==32)
elder<-subset(des,RIDAGEYR>=60);core<-subset(des,RIDAGEYR>=60 & common_scores==1);stopifnot(nrow(model.frame(core))==3336)
out<-list();k<-0
add<-function(model,term,label,group='pooled',scale='difference',family='exploratory',design=core) {
 b<-coef(model)[term]; se<-sqrt(diag(vcov(model)))[term]; df<-degf(design);p<-2*pt(abs(b/se),df=df,lower.tail=FALSE);ci<-b+c(-1,1)*qt(.975,df)*se
 if(scale %in% c('odds ratio','prevalence ratio')){b<-exp(b);ci<-exp(ci)}
 k<<-k+1;out[[k]]<<-data.frame(analysis=label,group=group,n=nrow(model.frame(model)),estimate=unname(b),se_link=unname(se),ci_low=unname(ci[1]),ci_high=unname(ci[2]),p_value=unname(p),df=df,scale=scale,family=family)
}
fit<-function(form,des,term,label,group='pooled',family='exploratory',fam=gaussian(),scale='difference'){
 m<-svyglm(form,design=des,family=fam,na.action=na.omit);add(m,term,label,group,scale,family,des);invisible(m)
}
basic<-c('fi_all','fi_no_sleep','fi_no_fatigue','fi_neither','fi_difference','fi_diff_nosleep','fi_diff_nofatigue')
for(y in basic){form<-as.formula(paste(y,'~sleep_told_doctor + age_c + sex_f + cycle_f'));fit(form,core,'sleep_told_doctor',y,family='primary')}

fullcase<-subset(core,!is.na(race_f)&!is.na(edu_f)&!is.na(INDFMPIR)&!is.na(ever_smoked))
fit(fi_neither~sleep_told_doctor+age_c+sex_f+cycle_f,fullcase,'sleep_told_doctor','fi_neither_fullcase_basic')
fit(fi_neither~sleep_told_doctor+age_c+sex_f+cycle_f+race_f+edu_f+INDFMPIR+ever_smoked,fullcase,'sleep_told_doctor','fi_neither_fullcase_adjusted')
for(cy in levels(d$cycle_f)){
 cd<-subset(des,cycle_f==cy);cdcore<-subset(cd,RIDAGEYR>=60&common_scores==1)
 fit(fi_neither~sleep_told_doctor+age_c+sex_f,cdcore,'sleep_told_doctor','fi_neither_cycle',cy)
}
fit(fi_neither~sleep_told_doctor*cycle_f+age_c+sex_f,core,'sleep_told_doctor:cycle_f2013-2014','sleep_by_cycle_interaction')
fit(fi_neither~sleep_told_doctor*sex_f+age_c+cycle_f,core,'sleep_told_doctor:sex_ffemale','sleep_by_sex_interaction')
for(y in c('fi_no_mood','fi_no_mood_no_care','fi_neither_no_care','fi_neither_no_labs','fi_neither_no_bmi','fi_neither_no_function_walk','fi_missing5_subset')){
 form<-as.formula(paste(y,'~sleep_told_doctor + age_c + sex_f + cycle_f'));fit(form,core,'sleep_told_doctor',y)
}
fit(fi_neither_complete~sleep_told_doctor+age_c+sex_f+cycle_f,elder,'sleep_told_doctor','fi_neither_complete_cases')
fit(fi_neither~sleep_told_doctor+age_c+sex_f+cycle_f,subset(core,RIDAGEYR>=65),'sleep_told_doctor','fi_neither_age65')
fit(fi_neither~sleep_recent_halfdays+age_c+sex_f+cycle_f,core,'sleep_recent_halfdays','recent_sleep_symptoms_ge_halfdays')
fit(fi_neither~sleep_recent_any+age_c+sex_f+cycle_f,core,'sleep_recent_any','recent_sleep_symptoms_any')
fit(fi_high~sleep_told_doctor+age_c+sex_f+cycle_f,core,'sleep_told_doctor','FI_above_0.21_OR',fam=quasibinomial(),scale='odds ratio')
fit(fi_high~sleep_told_doctor+age_c+sex_f+cycle_f,core,'sleep_told_doctor','FI_above_0.21_PR',fam=quasipoisson(link='log'),scale='prevalence ratio')
comp<-subset(core,!is.na(prior_physical4)&!is.na(fried_physical4));fit(fi_neither~sleep_told_doctor+age_c+sex_f+cycle_f,comp,'sleep_told_doctor','comparison_fi_neither');fit(prior_physical4~sleep_told_doctor+age_c+sex_f+cycle_f,comp,'sleep_told_doctor','comparison_prior_physical4');fit(fried_physical4~sleep_told_doctor+age_c+sex_f+cycle_f,comp,'sleep_told_doctor','comparison_fried_physical4')
comp<-update(comp,score_difference=fi_neither-prior_physical4);fit(score_difference~sleep_told_doctor+age_c+sex_f+cycle_f,comp,'sleep_told_doctor','paired_fi_minus_prior')
mort<-subset(elder,ELIGSTAT==1 & !is.na(MORTSTAT)&!is.na(fi_neither));cat('mortality n',nrow(model.frame(mort)),'events',sum(model.frame(mort)$MORTSTAT==1),'\n')
fit(MORTSTAT~fi_10+age_c+sex_f+cycle_f,mort,'fi_10','mortality_OR_per_0.1_FI',fam=quasibinomial(),scale='odds ratio')
fit(MORTSTAT~fi_all10+age_c+sex_f+cycle_f,mort,'fi_all10','mortality_OR_per_0.1_all_FI',fam=quasibinomial(),scale='odds ratio')
fit(MORTSTAT~fi_nomoodcare10+age_c+sex_f+cycle_f,mort,'fi_nomoodcare10','mortality_OR_per_0.1_no_mood_care',fam=quasibinomial(),scale='odds ratio')
fit(MORTSTAT~fi_10+age_c+sex_f+cycle_f,mort,'fi_10','mortality_PR_per_0.1_FI',fam=quasipoisson(link='log'),scale='prevalence ratio')
miss<-subset(elder,!is.na(sleep_told_doctor)&!is.na(race_f)&!is.na(edu_f)&!is.na(INDFMPIR)&!is.na(ever_smoked))
m_missing<-fit(I(1-FI_eligible)~sleep_told_doctor+age_c+sex_f+cycle_f+race_f+edu_f+INDFMPIR+ever_smoked,miss,'sleep_told_doctor','missing_FI_OR_sleep',fam=quasibinomial(),scale='odds ratio')
add(m_missing,'age_c','missing_FI_OR_per_year',scale='odds ratio',design=miss)
add(m_missing,'sex_ffemale','missing_FI_OR_female',scale='odds ratio',design=miss)
for(sex in c('male','female')){
 sub<-subset(core,sex_f==sex);fit(fi_neither~sleep_told_doctor+age_c+cycle_f,sub,'sleep_told_doctor','fi_neither_by_sex',sex)
}
# Refit propensity and outcome in every survey JKn replicate to include estimated response-weight uncertainty.
repdesign<-as.svrepdesign(des,type='JKn');repdomain<-subset(repdesign,RIDAGEYR>=60 & !is.na(sleep_told_doctor)&!is.na(race_f)&!is.na(edu_f)&!is.na(INDFMPIR)&!is.na(ever_smoked))
response_fun<-function(w,data){
 data$response<-as.integer(data$common_scores==1);w<-w/mean(w[w>0])
 propensity<-glm(response~sleep_told_doctor+age_c+sex_f+cycle_f+race_f+edu_f+INDFMPIR+ever_smoked,data=data,weights=w,family=quasibinomial())
 pr<-predict(propensity,type='response');keep<-data$response==1 & is.finite(pr) & pr>0 & w>0
 result<-lm(fi_neither~sleep_told_doctor+age_c+sex_f+cycle_f+race_f+edu_f+INDFMPIR+ever_smoked,data=data[keep,],weights=w[keep]/pr[keep])
 coef(result)[['sleep_told_doctor']]
}
ipw_data<-model.frame(repdomain);ipw_data$response<-as.integer(ipw_data$common_scores==1);ipw_w<-weights(repdomain,type='sampling');ipw_w<-ipw_w/mean(ipw_w[ipw_w>0]);ipw_model<-glm(response~sleep_told_doctor+age_c+sex_f+cycle_f+race_f+edu_f+INDFMPIR+ever_smoked,data=ipw_data,weights=ipw_w,family=quasibinomial());ipw_p<-predict(ipw_model,type='response');ipw_wobs<-ipw_w[ipw_data$response==1]/ipw_p[ipw_data$response==1];write.csv(data.frame(n_total=nrow(ipw_data),n_responders=sum(ipw_data$response),propensity_min=min(ipw_p),propensity_median=median(ipw_p),propensity_max=max(ipw_p),inverse_multiplier_max=max(1/ipw_p),kish_weight_ess=sum(ipw_wobs)^2/sum(ipw_wobs^2)), 'supplement_ipw_diagnostics.csv',row.names=FALSE)
ipw<-withReplicates(repdomain,response_fun);ipwb<-as.numeric(coef(ipw));ipwse<-as.numeric(SE(ipw));ipwdf<-degf(des);ipwp<-2*pt(abs(ipwb/ipwse),df=ipwdf,lower.tail=FALSE)
k<-k+1;out[[k]]<-data.frame(analysis='fi_neither_response_IPW_reestimated_JKn',group='pooled',n=sum(model.frame(repdomain)$common_scores==1),estimate=ipwb,se_link=ipwse,ci_low=ipwb-qt(.975,ipwdf)*ipwse,ci_high=ipwb+qt(.975,ipwdf)*ipwse,p_value=ipwp,df=ipwdf,scale='difference',family='exploratory')
means<-list();j<-0
getmean<-function(des,var,label){
 f<-as.formula(paste0('~',var));v<-svymean(f,des,na.rm=TRUE);b<-coef(v)[1];s<-SE(v)[1];df<-degf(des);j<<-j+1;means[[j]]<<-data.frame(label=label,n=sum(!is.na(model.frame(des)[[var]])),estimate=unname(b),se=unname(s),ci_low=unname(b-qt(.975,df)*s),ci_high=unname(b+qt(.975,df)*s),df=df)
}
for(v in c('RIDAGEYR','sleep_told_doctor','fi_all','fi_neither','fi_high','FI_eligible')){getmean(elder,v,paste0('elder_',v));getmean(core,v,paste0('core_',v))}
for(cy in levels(d$cycle_f)){
 ds<-subset(des,cycle_f==cy & RIDAGEYR>=60);getmean(ds,'fi_neither',paste0(cy,'_FI'));getmean(ds,'FI_eligible',paste0(cy,'_FI_eligible'))
}
for(gr in c('male','female'))getmean(subset(core,sex_f==gr),'fi_neither',paste0('core_FI_',gr))
for(gr in c(0,1))getmean(subset(core,sleep_told_doctor==gr),'fi_neither',paste0('core_FI_sleep_',gr))
for(gr in c(0,1))getmean(subset(mort,fi_high==gr),'MORTSTAT',paste0('mortality_FI_high_',gr))
tabbase<-update(core,female=as.numeric(RIAGENDR==2),race_mexican=as.numeric(RIDRETH3==1),race_other_hispanic=as.numeric(RIDRETH3==2),race_NHW=as.numeric(RIDRETH3==3),race_NHB=as.numeric(RIDRETH3==4),race_asian=as.numeric(RIDRETH3==6),race_other=as.numeric(RIDRETH3==7),education_below_highschool=as.numeric(education_clean<=2))
for(group in c('all','no_disclosure','disclosure')){
 ds<-if(group=='all')tabbase else if(group=='no_disclosure')subset(tabbase,sleep_told_doctor==0) else subset(tabbase,sleep_told_doctor==1)
 for(v in c('RIDAGEYR','female','race_mexican','race_other_hispanic','race_NHW','race_NHB','race_asian','race_other','education_below_highschool','INDFMPIR','ever_smoked','BMXBMI','fi_all','fi_neither'))getmean(ds,v,paste0('baseline_',group,'_',v))
}
out<-do.call(rbind,out);out$p_holm_family<-NA_real_
for(f in unique(out$family)){idx<-which(out$family==f);out$p_holm_family[idx]<-p.adjust(out$p_value[idx],method='holm')}
meansdf<-do.call(rbind,means);write.csv(out,'supplement_survey_models.csv',row.names=FALSE,na='');write.csv(meansdf,'supplement_weighted_descriptives.csv',row.names=FALSE,na='');write.csv(meansdf[grepl('^baseline_',meansdf$label),],'supplement_baseline_table.csv',row.names=FALSE,na='')
cat('R survey version',as.character(packageVersion('survey')),'\n');print(out[,c('analysis','group','n','estimate','ci_low','ci_high','p_value','df')],row.names=FALSE)

#!/usr/bin/env Rscript
.libPaths(c('/Volumes/Brain6NHANES/Rlib','/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
suppressPackageStartupMessages(library(survey))
out<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/phenotype_pilot/visit_code_audit/original'
old<-'/Volumes/Brain6Sim/nhanes-public-2011-14-20261028';new<-'/Volumes/Brain6NHANES/source'
readfile<-function(name,s){path<-file.path(if(name %in% c('DEMO','DPQ','PFQ','MGX'))old else new,paste0(name,'_',s,'.XPT'));foreign::read.xport(path)}
code<-function(x,m){y<-rep(NA_real_,length(x));for(k in names(m))y[!is.na(x)&x==as.numeric(k)]<-m[[k]];y}
yesno<-function(x)code(x,c('1'=1,'2'=0));difficulty<-function(x)code(x,c('1'=0,'2'=1,'3'=1,'4'=1))
missingmean<-function(x,min_n){n<-rowSums(!is.na(x));v<-rowMeans(x,na.rm=TRUE);v[n<min_n]<-NA_real_;v}
data_all<-list();fi_items<-NULL
for(s in c('G','H')){
 d<-readfile('DEMO',s);stopifnot(!anyDuplicated(d$SEQN))
 for(nm in c('DPQ','PFQ','MGX','MCQ','BPQ','DIQ','KIQ_U','HUQ','BMX','GHB','CBC','PAQ','WHQ')){
  z<-readfile(nm,s);stopifnot(!anyDuplicated(z$SEQN));d<-merge(d,z,by='SEQN',all.x=TRUE,sort=FALSE)
 }
 rx<-readfile('RXQ_RX',s);groups<-split(rx,rx$SEQN)
 meds<-do.call(rbind,lapply(groups,function(g){u<-unique(g$RXDUSE[!is.na(g$RXDUSE)]);v<-unique(g$RXDCOUNT[!is.na(g$RXDCOUNT)]);
  stopifnot(length(u)<=1,length(v)<=1);data.frame(SEQN=g$SEQN[1],use=if(length(u))u else NA_real_,count=if(length(v))v else NA_real_)}))
 d<-merge(d,meds,by='SEQN',all.x=TRUE,sort=FALSE);stopifnot(!anyDuplicated(d$SEQN))
 d<-d[!is.na(d$WTMEC2YR) & d$WTMEC2YR>0 & !is.na(d$SDMVSTRA) & !is.na(d$SDMVPSU),]
 n<-nrow(d);bmi<-d$BMXBMI;sex<-d$RIAGENDR
 cpt<-list();cpt$bmi<-ifelse(is.na(bmi),NA,ifelse(bmi<18.5 | bmi>=30,1,ifelse(bmi<=25,0,.5)))
 gripcols<-paste0('MGXH',rep(1:2,each=3),'T',rep(1:3,2));trials<-as.matrix(d[,gripcols]);trials[!is.na(trials)&trials<=0]<-NA
 grip<-apply(trials,1,function(x)if(all(is.na(x))) NA_real_ else max(x,na.rm=TRUE))
 lim<-ifelse(sex==1,ifelse(bmi<=24,29,ifelse(bmi<=28,30,32)),ifelse(sex==2,ifelse(bmi<=23,17,ifelse(bmi<=26,17.3,ifelse(bmi<=29,18,21))),NA))
 cpt$grip<-ifelse(is.na(grip)|is.na(lim),NA,as.numeric(grip<=lim))
 for(v in c('MCQ160K','MCQ160M','MCQ160A','MCQ160B','MCQ220','MCQ160C','MCQ160D','KIQ022','MCQ160E','BPQ020','MCQ160F'))cpt[[v]]<-yesno(d[[v]])
 leak<-d[,c('KIQ042','KIQ044','KIQ046')];anyyes<-rowSums(leak==1,na.rm=TRUE)>0;allno<-rowSums(leak==2,na.rm=TRUE)==3
 cpt$urinary<-ifelse(anyyes,1,ifelse(allno,0,NA));cpt$diabetes<-code(d$DIQ010,c('1'=1,'2'=0,'3'=.5));cpt$confusion<-yesno(d$PFQ057)
 cpt$medications<-ifelse(d$use==2,0,ifelse(d$use==1 & !is.na(d$count),ifelse(d$count>=5,1,.5),NA))
 cpt$hospital<-yesno(d$HUQ071);cpt$health_change<-code(d$HUQ020,c('1'=0,'2'=1,'3'=0))
 cpt$health<-code(d$HUQ010,c('1'=0,'2'=0,'3'=0,'4'=1,'5'=1))
 visits<-if(s=='G')d$HUQ050 else d$HUQ051
 cpt$visits<-if(s=='G')code(visits,c('0'=0,'1'=.5,'2'=.5,'4'=1,'5'=1)) else code(visits,c('0'=0,'1'=.5,'2'=.5,'3'=.5,'4'=1,'5'=1,'6'=1,'7'=1))
 pfq<-c('A','B','C','D','E','F','R','G','M','I','L','K','P','J','T')
 for(k in pfq)cpt[[paste0('PFQ',k)]]<-difficulty(d[[paste0('PFQ061',k)]])
 labs<-list(LBXGH=c(-Inf,5.7),LBXNEPCT=c(40,80),LBXLYPCT=c(20,40),LBXRDW=c(11.6,14.6))
 for(k in names(labs)){x<-d[[k]];bounds<-labs[[k]];cpt[[k]]<-ifelse(is.na(x),NA,as.numeric(x<bounds[1] | x>bounds[2]))}
 for(k in c('LBXHGB','LBXRBCSI')){x<-d[[k]];lo<-ifelse(sex==1,if(k=='LBXHGB')13.5 else 4.7,ifelse(sex==2,if(k=='LBXHGB')12 else 4.2,NA));hi<-ifelse(sex==1,if(k=='LBXHGB')18 else 6.1,ifelse(sex==2,if(k=='LBXHGB')16 else 5.4,NA));cpt[[k]]<-ifelse(is.na(x)|is.na(lo),NA,as.numeric(x<lo | x>hi))}
 for(k in c('DPQ030','DPQ070','DPQ010','DPQ040','DPQ020','DPQ050','DPQ060'))cpt[[k]]<-code(d[[k]],c('0'=0,'1'=.33,'2'=.66,'3'=1))
 stopifnot(length(cpt)==49,!anyDuplicated(names(cpt)));mat<-do.call(cbind,cpt);fi_items<-colnames(mat)
 stopifnot(all(is.na(mat) | mat %in% c(0,.33,.5,.66,1)))
 # Five GPAQ domains; no activity means zero; invalid/unknown days or duration remain missing.
 activity<-list(c('PAQ605','PAQ610','PAD615',8),c('PAQ620','PAQ625','PAD630',4),c('PAQ635','PAQ640','PAD645',4),c('PAQ650','PAQ655','PAD660',8),c('PAQ665','PAQ670','PAD675',4))
 domain<-lapply(activity,function(a){yes<-d[[a[1]]];days<-d[[a[2]]];mins<-d[[a[3]]];v<-rep(NA_real_,n);v[!is.na(yes)&yes==2]<-0;ok<-!is.na(yes)&yes==1&!is.na(days)&days>=1&days<=7&!is.na(mins)&mins>=1&mins<=1440;v[ok]<-days[ok]*mins[ok]*as.numeric(a[4]);v})
 met<-rowSums(do.call(cbind,domain),na.rm=FALSE);lowactivity<-ifelse(is.na(met),NA,as.numeric(met<600))
 validwt<-!is.na(d$WHD020)&!is.na(d$WHD050)&d$WHD020>0&d$WHD020<7777&d$WHD050>0&d$WHD050<7777
 loss<-d$WHD050-d$WHD020;loss5<-validwt & loss>=.05*d$WHD050;intent<-d$WHQ060
 lossflag<-ifelse(!validwt,NA,ifelse(!loss5,0,ifelse(intent==2,1,ifelse(intent==1,0,NA))))
 bmiflag<-ifelse(is.na(bmi),NA,as.numeric(bmi<=22.5))
 wloss<-ifelse((!is.na(bmiflag)&bmiflag==1)|(!is.na(lossflag)&lossflag==1),1,
               ifelse(!is.na(bmiflag)&bmiflag==0 & !is.na(lossflag)&lossflag==0,0,NA))
 for(mode in c('primary','mobility_proxy','visits_low','visits_high')){
  mm<-mat
  if(mode=='mobility_proxy')for(k in c('PFQB','PFQC')){sel<-is.na(mm[,k])&!is.na(d$PFQ054)&d$PFQ054==1;mm[sel,k]<-1}
  if(s=='G' && mode %in% c('visits_low','visits_high'))mm[!is.na(visits)&visits==3,'visits']<-if(mode=='visits_low').5 else 1
  full<-missingmean(mm,40);fi47<-missingmean(mm[,!colnames(mm)%in%c('DPQ030','DPQ040'),drop=FALSE],0)
  fi48<-missingmean(mm[,colnames(mm)!='DPQ030',drop=FALSE],0)
  fi_noshared<-missingmean(mm[,!colnames(mm)%in%c('PFQB','PFQE'),drop=FALSE],0)
  fi47_noshared<-missingmean(mm[,!colnames(mm)%in%c('DPQ030','DPQ040','PFQB','PFQE'),drop=FALSE],0)
  fi43<-missingmean(mm[,!colnames(mm)%in%c('DPQ030','DPQ040','PFQB','PFQE','bmi','grip'),drop=FALSE],0)
  z<-!is.na(mm[,'DPQ030']) & !is.na(mm[,'DPQ040']) & rowSums(!is.na(mm))>=40
  check47<-(rowSums(mm,na.rm=TRUE)-mm[,'DPQ030']-mm[,'DPQ040'])/(rowSums(!is.na(mm))-2)
  stopifnot(all(abs(fi47[z]-check47[z])<1e-12))
  ffp<-cbind(wloss,difficulty(d$PFQ061H),difficulty(d$PFQ061E),mm[,'PFQB'],lowactivity)
  ffs<-rowSums(ffp,na.rm=FALSE);ffs4<-rowSums(ffp[,-4,drop=FALSE],na.rm=FALSE)
  ffs3<-rowSums(ffp[,c(1,2,5),drop=FALSE],na.rm=FALSE);ffs2<-rowSums(ffp[,c(2,5),drop=FALSE],na.rm=FALSE)
  stopifnot(all(is.na(ffs) | ffs %in% 0:5),all(is.na(ffs4) | ffs4 %in% 0:4))
  data_all[[paste(s,mode)]]<-data.frame(cycle=s,mode=mode,SEQN=d$SEQN,STRATUM=paste(s,d$SDMVSTRA,sep='_'),
   SDMVPSU=d$SDMVPSU,weight=d$WTMEC2YR/2,age=d$RIDAGEYR,fi=full,fi48=fi48,fi47=fi47,fi_noshared=fi_noshared,fi47_noshared=fi47_noshared,fi43=fi43,ffs=ffs,ffs4=ffs4,ffs3=ffs3,ffs2=ffs2,
   sleep=mm[,'DPQ030'],fatigue=mm[,'DPQ040'])
 }
}
writeLines(fi_items,file.path(out,'fi_components.txt'))
correlate<-function(design,xname,yname){
 f<-as.formula(paste0('~',xname,'+',yname,'+I(',xname,'^2)+I(',yname,'^2)+I(',xname,'*',yname,')'))
 m<-svymean(f,design)
 # Rename both coefficient vector and covariance dimensions for a stable contrast.
 names(m)<-c('x','y','xx','yy','xy');dimnames(attr(m,'var'))<-list(names(m),names(m))
 c<-svycontrast(m,quote((xy-x*y)/sqrt((xx-x*x)*(yy-y*y))))
 estimate<-unname(coef(c));se<-unname(SE(c));critical<-qt(.975,df=degf(design))
 z<-as.numeric(m);fcheck<-function(t)(t[5]-t[1]*t[2])/sqrt((t[3]-t[1]^2)*(t[4]-t[2]^2))
 grad<-vapply(seq_along(z),function(j){h<-rep(0,length(z));h[j]<-1e-6;(fcheck(z+h)-fcheck(z-h))/(2e-6)},numeric(1))
 checkse<-sqrt(drop(t(grad)%*%vcov(m)%*%grad));stopifnot(abs(checkse-se)<1e-7)
 c(estimate=estimate,se=se,low=estimate-critical*se,high=estimate+critical*se)
}
results<-list();counts<-list();paired<-list()
for(mode in c('primary','mobility_proxy','visits_low','visits_high')){
 a<-do.call(rbind,data_all[grep(paste0(' ',mode,'$'),names(data_all))]);stopifnot(!anyDuplicated(a$SEQN))
 a$eligible<-!is.na(a$age)&a$age>=60 & complete.cases(a[,c('fi','fi48','fi47','ffs','ffs4','sleep','fatigue')])
 design<-svydesign(id=~SDMVPSU,strata=~STRATUM,weights=~weight,nest=TRUE,data=a)
 sub<-subset(design,eligible);counts[[length(counts)+1]]<-data.frame(mode=mode,design_n=nrow(a),eligible_n=sum(a$eligible),df=degf(sub))
 for(pair in list(c('fi','ffs'),c('fi','sleep'),c('fi47','sleep'),c('fi47','ffs4'))){
  v<-correlate(sub,pair[1],pair[2]);q<-a[a$eligible,];x<-q[[pair[1]]];y<-q[[pair[2]]];w<-q$weight
  mx<-sum(w*x)/sum(w);my<-sum(w*y)/sum(w)
  direct<-sum(w*(x-mx)*(y-my))/sqrt(sum(w*(x-mx)^2)*sum(w*(y-my)^2))
  stopifnot(abs(direct-v['estimate'])<1e-12)
  results[[length(results)+1]]<-data.frame(mode=mode,x=pair[1],y=pair[2],t(v))
 }
 m<-svymean(~fi+fi47+sleep+I(fi^2)+I(fi47^2)+I(sleep^2)+I(fi*sleep)+I(fi47*sleep),sub)
 names(m)<-c('a','b','s','aa','bb','ss','u','v');dimnames(attr(m,'var'))<-list(names(m),names(m))
 contrast<-svycontrast(m,quote((u-a*s)/sqrt((aa-a*a)*(ss-s*s))-(v-b*s)/sqrt((bb-b*b)*(ss-s*s))))
 est<-unname(coef(contrast));se<-unname(SE(contrast));crit<-qt(.975,df=degf(sub))
 z<-as.numeric(m);fdiff<-function(t)(t[7]-t[1]*t[3])/sqrt((t[4]-t[1]^2)*(t[6]-t[3]^2))-(t[8]-t[2]*t[3])/sqrt((t[5]-t[2]^2)*(t[6]-t[3]^2))
 grad<-vapply(seq_along(z),function(j){h<-rep(0,length(z));h[j]<-1e-6;(fdiff(z+h)-fdiff(z-h))/(2e-6)},numeric(1))
 stopifnot(abs(sqrt(drop(t(grad)%*%vcov(m)%*%grad))-se)<1e-7)
 paired[[length(paired)+1]]<-data.frame(mode=mode,comparison='corr(fi,sleep)-corr(fi47,sleep)',estimate=est,se=se,low=est-crit*se,high=est+crit*se)
}
write.table(do.call(rbind,counts),file.path(out,'counts.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,results),file.path(out,'correlations.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,paired),file.path(out,'paired_contrasts.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('Completed exploratory survey-weighted scoring; no participant-level outputs saved.\n');print(sessionInfo())

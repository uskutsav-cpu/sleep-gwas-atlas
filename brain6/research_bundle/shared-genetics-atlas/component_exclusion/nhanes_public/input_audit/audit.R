#!/usr/bin/env Rscript
source_dir<-'/Volumes/Brain6Sim/nhanes-public-2011-14-20261028'
out<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/input_audit'
cols<-list(DEMO=c('SEQN','RIDAGEYR','WTMEC2YR','SDMVPSU','SDMVSTRA'),
 DPQ=c('SEQN','DPQ030','DPQ040'),PFQ=c('SEQN','PFQ054','PFQ061B','PFQ061E'),
 MGX=c('SEQN','MGXH1T1','MGXH2T1','MGXH1T2','MGXH2T2','MGXH1T3','MGXH2T3'))
rows<-list();inventory<-list()
for(suffix in c('G','H')){
 data<-list()
 for(name in names(cols)){
  x<-foreign::read.xport(file.path(source_dir,paste0(name,'_',suffix,'.XPT')))
  stopifnot(all(cols[[name]] %in% names(x)),!anyDuplicated(x$SEQN))
  inventory[[length(inventory)+1]]<-data.frame(cycle=suffix,file=name,rows=nrow(x),columns=ncol(x),required_present=TRUE)
  data[[name]]<-x[,cols[[name]],drop=FALSE]
 }
 d<-data$DEMO
 for(name in c('DPQ','PFQ','MGX'))d<-merge(d,data[[name]],by='SEQN',all.x=TRUE,sort=FALSE)
 stopifnot(nrow(d)==nrow(data$DEMO),!anyDuplicated(d$SEQN))
 a<-d[d$RIDAGEYR>=60,,drop=FALSE];valid<-function(x,values)!is.na(x)&x %in% values
 fields<-c('DPQ030','DPQ040','PFQ061B','PFQ061E')
 for(field in fields){
  ok<-valid(a[[field]],if(startsWith(field,'DPQ'))0:3 else 1:4)
  rows[[length(rows)+1]]<-data.frame(cycle=suffix,metric=paste0(field,'_valid'),n=sum(ok),denominator=nrow(a))
 }
 all4<-valid(a$DPQ030,0:3)&valid(a$DPQ040,0:3)&valid(a$PFQ061B,1:4)&valid(a$PFQ061E,1:4)
 grip<-rowSums(!is.na(a[,cols$MGX[-1],drop=FALSE]) & a[,cols$MGX[-1],drop=FALSE]>0)>0
 checks<-list(joint_four_valid=all4,at_least_one_positive_grip=grip,joint_four_grip=all4 & grip,
              eligible_design=(!is.na(a$WTMEC2YR)&a$WTMEC2YR>0 & !is.na(a$SDMVPSU)&!is.na(a$SDMVSTRA)),
              mobility_equipment_yes=valid(a$PFQ054,1),
              mobility_equipment_yes_walking_missing=valid(a$PFQ054,1)&is.na(a$PFQ061B),
              mobility_equipment_no_walking_missing=valid(a$PFQ054,2)&is.na(a$PFQ061B))
 for(nm in names(checks))rows[[length(rows)+1]]<-data.frame(cycle=suffix,metric=nm,n=sum(checks[[nm]]),denominator=nrow(a))
}
write.table(do.call(rbind,rows),file.path(out,'counts.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,inventory),file.path(out,'files.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS: required columns present; unique participant keys; no person-level outputs written.\n')

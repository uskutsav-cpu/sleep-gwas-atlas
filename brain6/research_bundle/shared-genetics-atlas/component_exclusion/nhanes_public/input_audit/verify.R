#!/usr/bin/env Rscript
root<-'/Volumes/Brain6Sim/nhanes-public-2011-14-20261028'
out<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/input_audit'
reported<-read.delim(file.path(out,'counts.tsv'))
for(s in c('G','H')){
 d<-foreign::read.xport(file.path(root,paste0('DEMO_',s,'.XPT')))
 p<-foreign::read.xport(file.path(root,paste0('PFQ_',s,'.XPT')))
 q<-foreign::read.xport(file.path(root,paste0('DPQ_',s,'.XPT')))
 i<-which(d$RIDAGEYR>=60);ids<-d$SEQN[i];p<-p[match(ids,p$SEQN),];q<-q[match(ids,q$SEQN),]
 vals<-c(DPQ030_valid=sum(q$DPQ030 %in% 0:3),DPQ040_valid=sum(q$DPQ040 %in% 0:3),
   PFQ061B_valid=sum(p$PFQ061B %in% 1:4),PFQ061E_valid=sum(p$PFQ061E %in% 1:4),
   joint_four_valid=sum(q$DPQ030 %in% 0:3 & q$DPQ040 %in% 0:3 & p$PFQ061B %in% 1:4 & p$PFQ061E %in% 1:4),
   mobility_equipment_yes=sum(p$PFQ054 %in% 1),
   mobility_equipment_yes_walking_missing=sum(p$PFQ054 %in% 1 & is.na(p$PFQ061B)),
   mobility_equipment_no_walking_missing=sum(p$PFQ054 %in% 2 & is.na(p$PFQ061B)))
 for(k in names(vals)){r<-reported[reported$cycle==s & reported$metric==k,];stopifnot(nrow(r)==1,r$n==vals[k],r$denominator==length(ids))}
}
cat('PASS: sleep, fatigue, physical-item and skip counts verified with match-based joins.\n')

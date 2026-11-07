#!/usr/bin/env Rscript
out<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/item_mapping'
src<-'/Volumes/Brain6NHANES/source';demo<-'/Volumes/Brain6Sim/nhanes-public-2011-14-20261028'
reported<-read.delim(file.path(out,'ambiguity_counts.tsv'))
for(s in c('G','H')){
 d<-foreign::read.xport(file.path(demo,paste0('DEMO_',s,'.XPT')))
 h<-foreign::read.xport(file.path(src,paste0('HUQ_',s,'.XPT')))
 w<-foreign::read.xport(file.path(src,paste0('WHQ_',s,'.XPT')))
 b<-foreign::read.xport(file.path(src,paste0('BMX_',s,'.XPT')))
 v<-if(s=='G')'HUQ050' else 'HUQ051'
 a<-merge(merge(merge(d[,c('SEQN','RIDAGEYR')],h[,c('SEQN',v)],by='SEQN',all.x=TRUE),
  w[,c('SEQN','WHD020','WHD050','WHQ060')],by='SEQN',all.x=TRUE),b[,c('SEQN','BMXBMI')],by='SEQN',all.x=TRUE)
 a<-subset(a,RIDAGEYR>=60); valid<-with(a,!is.na(WHD020)&!is.na(WHD050)&WHD020>0&WHD020<7777&WHD050>0&WHD050<7777)
 candidate<-subset(a,valid & !is.na(BMXBMI) & BMXBMI>22.5 & WHD020 <= .95*WHD050 & WHD050-WHD020<10 & is.na(WHQ060))
 checks<-c(older_adults=nrow(a),visit_category_3=sum(a[[v]] %in% 3),valid_self_report_weights=sum(valid),weight_loss_criterion_intent_unknown=nrow(candidate))
 for(k in names(checks)){r<-reported[reported$cycle==s & reported$metric==k,];stopifnot(nrow(r)==1,r$n==checks[k])}
}
cat('PASS: matched older-adult visit and weight ambiguity counts independently checked.\n')

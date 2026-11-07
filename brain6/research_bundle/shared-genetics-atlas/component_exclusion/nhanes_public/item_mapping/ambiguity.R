#!/usr/bin/env Rscript
src<-'/Volumes/Brain6NHANES/source';demo<-'/Volumes/Brain6Sim/nhanes-public-2011-14-20261028'
out<-'outputs/shared-genetics-atlas/component_exclusion/nhanes_public/item_mapping';result<-list()
for(s in c('G','H')){
 d<-foreign::read.xport(file.path(demo,paste0('DEMO_',s,'.XPT')))[,c('SEQN','RIDAGEYR')]
 h<-foreign::read.xport(file.path(src,paste0('HUQ_',s,'.XPT')))
 w<-foreign::read.xport(file.path(src,paste0('WHQ_',s,'.XPT')))
 b<-foreign::read.xport(file.path(src,paste0('BMX_',s,'.XPT')))
 stopifnot(all(vapply(list(d,h,w,b),function(x)!anyDuplicated(x$SEQN),logical(1))))
 d<-d[d$RIDAGEYR>=60,];h<-h[match(d$SEQN,h$SEQN),];w<-w[match(d$SEQN,w$SEQN),];b<-b[match(d$SEQN,b$SEQN),]
 visits<-if(s=='G')h$HUQ050 else h$HUQ051
 validwt<-!is.na(w$WHD020)&!is.na(w$WHD050)&w$WHD020>0&w$WHD020<7777&w$WHD050>0&w$WHD050<7777
 loss<-w$WHD050-w$WHD020
 risk<-validwt & !is.na(b$BMXBMI) & b$BMXBMI>22.5 & loss>=.05*w$WHD050 & loss<10 & is.na(w$WHQ060)
 metrics<-c(older_adults=nrow(d),visit_category_3=sum(visits %in% 3),
  valid_self_report_weights=sum(validwt),weight_loss_criterion_intent_unknown=sum(risk,na.rm=TRUE))
 result[[length(result)+1]]<-data.frame(cycle=s,metric=names(metrics),n=unname(metrics))
}
write.table(do.call(rbind,result),file.path(out,'ambiguity_counts.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS: post-hoc unweighted count diagnostic; no person-level rows exported.\n')

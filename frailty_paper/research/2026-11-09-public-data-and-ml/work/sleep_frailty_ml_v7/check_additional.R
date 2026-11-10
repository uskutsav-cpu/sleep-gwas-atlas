# Independent base-R calculations for the added source. No model fitting.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
d <- read.csv(args[1],check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(d)==1600,!anyDuplicated(paste(d$source_ID,d$condition)))
stopifnot(sum(d$condition=='baseline')==547,sum(d$condition=='IL6')==517,sum(d$condition=='IFNa')==536)
out <- list()
for (condition in c('baseline','IL6','IFNa')) {
  x <- d[d$condition==condition,]
  y <- x$absolute_log2_effect
  for (arm in c('sequence_only','AVI_added')) {
    p <- x[[paste0(arm,'_baseline_context_projection')]]
    lo <- x[[paste0(arm,'_lower90')]]
    hi <- x[[paste0(arm,'_upper90')]]
    keep <- x[[paste0(arm,'_authorized_for_condition')]]==1
    groups <- split(seq_len(nrow(x)),x$local_group)
    groups <- groups[vapply(groups,function(z)length(z)>=5 && diff(range(y[z]))>0,TRUE)]
    local <- vapply(groups,function(z)if(diff(range(p[z]))>0)cor(y[z],p[z],method='spearman')else 0,0.)
    if(condition!='baseline')stopifnot(!any(keep))
    out[[length(out)+1]] <- data.frame(arm=arm,condition=condition,n=nrow(x),MAE=mean(abs(y-p)),RMSE=sqrt(mean((y-p)^2)),baseline_MAE=mean(abs(y-x$baseline_cell_median)),global_Spearman=if(diff(range(p))>0)cor(y,p,method='spearman')else 0,local_mean_Spearman=if(length(local))mean(local)else 0,local_groups=length(local),interval_coverage=mean(y>=lo&y<=hi),median_interval_width=median(hi-lo),ranked_fraction=mean(keep),ranked_observations=sum(keep))
  }
}
write.csv(do.call(rbind,out),args[2],row.names=FALSE)
cat('Independent additional-source arithmetic written; unsupported stimulation has no authorized answers\n')

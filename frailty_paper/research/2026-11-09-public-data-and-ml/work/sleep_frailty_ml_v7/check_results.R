# Independent base-R arithmetic from locked predictions and source-matched effects.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
d <- read.csv(args[1], check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(nrow(d)==9602, !anyDuplicated(paste(d$design_index,d$cell)))
out <- list()
for (arm in c('sequence_only','AVI_added')) {
  for (cell in sort(unique(d$cell))) {
    x <- d[d$cell==cell,]
    y <- x$measured_absolute_log2_effect
    p <- x[[paste0(arm,'_predicted_absolute_log2_effect')]]
    lo <- x[[paste0(arm,'_lower90')]]
    hi <- x[[paste0(arm,'_upper90')]]
    keep <- x[[paste0(arm,'_rankable')]]==1
    local <- split(seq_len(nrow(x)), x$local_group)
    local <- local[vapply(local,function(z)length(z)>=5 && diff(range(y[z]))>0,TRUE)]
    lr <- vapply(local,function(z)if(diff(range(p[z]))>0)cor(y[z],p[z],method='spearman') else 0,0.)
    out[[length(out)+1]] <- data.frame(arm=arm,cell=cell,n=nrow(x),MAE=mean(abs(y-p)),RMSE=sqrt(mean((y-p)^2)),baseline_MAE=mean(abs(y-x$baseline_cell_median)),global_Spearman=cor(y,p,method='spearman'),local_mean_Spearman=if(length(lr))mean(lr)else 0,local_groups=length(lr),interval_coverage=mean(y>=lo & y<=hi),median_interval_width=median(hi-lo),ranked_fraction=mean(keep),ranked_observations=sum(keep))
  }
}
write.csv(do.call(rbind,out),args[2],row.names=FALSE)
cat('Independent R numerical checks written:',nrow(d),'source-matched observations\n')

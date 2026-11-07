#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1)stop('Usage: Rscript plot_expanded.R RESULTS_DIR')
d<-read.delim(file.path(args[1], 'expanded_variant_evidence.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
if(nrow(d)!=425 || any(d$tissue_prediction_status!='PREDICTED'))stop('Plot requires all 425 completed variants')
ids<-unique(d$candidate_locus_id);stopifnot(length(ids)==6)
svg(file.path(args[1],'dnase_six_regions.svg'),width=12,height=8,pointsize=10)
par(mfrow=c(3,2),mar=c(4.2,4.5,2.5,1),las=1)
for (cid in ids) {
 x<-d[d$candidate_locus_id==cid,];ch<-unique(x$chromosome);stopifnot(length(ch)==1)
 xx<-as.numeric(x$pos_grch38_1based)/1e6; yy<-as.numeric(x$dnase_log2_ratio_501bp)
 lo<-min(c(0,yy));hi<-max(c(0,yy));pad<-max((hi-lo)*.12,.01)
 plot(xx,yy,pch=19,cex=.6,col=adjustcolor('#17618f',alpha.f=.75),xlab=paste(ch,'position, GRCh38 (Mb)'),ylab='Predicted DNase log2 ratio',main=paste0(sub('placo_insomnia__adhd_', '', cid),' (GRCh37), n = ',nrow(x)),ylim=c(lo-pad,hi+pad),bty='l')
 abline(h=0,col='#b87440',lty=2)
}
dev.off()

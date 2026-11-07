#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1)stop('Usage: Rscript plot_jbrowse.R OUTPUT_DIR')
d<-read.delim(file.path(args[1],'track_window.tsv'),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(d)==8004)
svg(file.path(args[1],'predicted_tracks.svg'),width=10,height=8,pointsize=10)
par(mfrow=c(4,1),mar=c(3,6.5,2.2,1),las=1)
for (j in list(c('dnase','0'),c('rna_seq','0'),c('rna_seq','1'),c('rna_seq','2'))) {
 x<-d[d$output_type==j[1] & d$track_ordinal==as.integer(j[2]),]
 stopifnot(nrow(x)==2001)
 name<-paste(x$track_name[1],x$life_stage[1],x$strand[1],sep=' | ')
 ylim<-range(c(x$reference,x$alternate));if(diff(ylim)==0)ylim<-ylim+c(-.001,.001)
 plot(x$offset_from_variant_bp,x$reference,type='l',col='#195e9e',ylim=ylim,xlab='Offset from rs2431108 (bp)',ylab='',main=name, yaxt='n')
 ticks<-pretty(ylim,n=3); ticks<-ticks[ticks>=min(ylim) & ticks<=max(ylim)]
 axis(2,at=ticks,labels=formatC(ticks,format='e',digits=1),las=1,cex.axis=.8)
 lines(x$offset_from_variant_bp,x$alternate,col='#d35f2e',lty=2)
 abline(v=0,col='gray40',lty=3)
 if(j[1]=='dnase')legend('topright',c('Reference T','Alternate C'),col=c('#195e9e','#d35f2e'),lty=c(1,2),bty='n',horiz=TRUE)
}
dev.off()

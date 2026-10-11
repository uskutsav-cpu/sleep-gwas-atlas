args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Usage: Rscript plot_transport.R source_table.tsv output_prefix")
d <- read.delim(args[1], stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(nrow(d)==3, all(d$ci_lower_95 < d$rg), all(d$ci_upper_95 > d$rg))
labels <- c("UKB complaints / Pan-UKB pain", "UKB complaints / MVP pain", "FinnGen registry insomnia / MVP pain")
cols <- c("#536779", "#536779", "#AC6200")
draw <- function() {
  par(mar=c(7.5,17,4,7), family="sans", xaxs="i", yaxs="i")
  plot(NA,xlim=c(0,0.7),ylim=c(0.5,3.5),xlab="Marginal genetic correlation",ylab="",yaxt="n",bty="n",cex.lab=1.05)
  abline(v=seq(0,0.7,0.1),col="#E8EBEE",lwd=0.8)
  axis(2,at=3:1,labels=labels,las=1,tick=FALSE,cex.axis=0.92)
  for (i in 1:3) {
    segments(d$ci_lower_95[i],4-i,d$ci_upper_95[i],4-i,col=cols[i],lwd=2.5)
    points(d$rg[i],4-i,pch=16,cex=1.35,col=cols[i])
    text(0.72,4-i,sprintf("%.3f\n[%.3f, %.3f]",d$rg[i],d$ci_lower_95[i],d$ci_upper_95[i]),adj=c(0,0.5),xpd=NA,cex=0.82,col=cols[i])
  }
  title("Source-defined insomnia and abdominal pain",adj=0,font.main=2,cex.main=1.18)
  mtext("Points and marginal 95% intervals; new FinnGen estimate is exploratory",side=3,line=0.4,adj=0,cex=0.86)
  mtext("Two MVP rows share an outcome GWAS. No calibrated between-source difference test.",side=1,line=4.6,adj=0,cex=0.79)
  mtext("Clinical/source constructs differ; CI95/log-Wald and constant effective-N assumptions remain.",side=1,line=5.8,adj=0,cex=0.79)
}
svg(paste0(args[2],".svg"),width=12,height=6.2);draw();dev.off()
png(paste0(args[2],".png"),width=2400,height=1240,res=200,type="cairo");draw();dev.off()

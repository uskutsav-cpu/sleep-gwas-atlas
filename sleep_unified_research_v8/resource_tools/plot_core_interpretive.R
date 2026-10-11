args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Usage: Rscript plot_core_interpretive.R core_rg.tsv output_prefix")
d <- read.delim(args[1], stringsAsFactors=FALSE, check.names=FALSE)
sleep <- unique(d$sleep_trait); outcome <- unique(d$outcome_trait)
stopifnot(nrow(d)==396, length(sleep)==12, length(outcome)==33)
held <- d$directional_interpretation_status=="WITHHOLD_PENDING_SOURCE_EFFECT_DIRECTION_RESOLUTION"
qc <- d$pairwise_h2_z_ge_4_both=="True" & d$pairwise_h2_intercept_le_1_2_both=="True"
stopifnot(sum(held)==36, all(is.na(d$rg[held])))
eligible <- !held & qc
limit <- max(abs(d$rg[eligible]))
limit <- ceiling(limit*10)/10
pal <- colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(101)
draw <- function() {
  par(mar=c(11.5,12,4,1), family="sans", xaxs="i", yaxs="i")
  plot(NA,xlim=c(0.5,33.5),ylim=c(0.5,12.5),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
  for(i in seq_len(nrow(d))) {
    x <- match(d$outcome_trait[i],outcome); y <- 13-match(d$sleep_trait[i],sleep)
    if(held[i]) {
      rect(x-.5,y-.5,x+.5,y+.5,col="#B8BDC5",border="white")
      text(x,y,"H",cex=.55,col="#252525")
    } else if(!qc[i]) {
      rect(x-.5,y-.5,x+.5,y+.5,col="white",border="#E2E2E2")
      text(x,y,"x",cex=.65,col="#555555")
    } else {
      k <- round((d$rg[i]+limit)/(2*limit)*100)+1
      rect(x-.5,y-.5,x+.5,y+.5,col=pal[k],border="white")
    }
  }
  axis(1,at=1:33,labels=gsub("_"," ",outcome),las=2,tick=FALSE,cex.axis=.76)
  axis(2,at=12:1,labels=gsub("_"," ",sleep),las=1,tick=FALSE,cex.axis=.8)
  title("Core source-defined genetic correlations: interpretation view",adj=0,cex.main=1.14)
  mtext(sprintf("396 historical records; 36 source-direction holds; %d further records hidden by pairwise QC",sum(!held & !qc)),side=3,line=.4,adj=0,cex=.83)
  legend("bottomleft",inset=c(0,-.27),xpd=NA,horiz=TRUE,bty="n",cex=.72,
         legend=c(sprintf("Negative rg (scale minimum %.1f)",-limit),"Near zero",sprintf("Positive rg (scale maximum %.1f)",limit),"H: source interpretation held","x: pairwise QC condition failed"),
         fill=c(pal[1],pal[51],pal[101],"#B8BDC5","white"),border=c(NA,NA,NA,NA,"#888888"))
  mtext("Colors show point estimates only; uncertainty is supplied in the release table. No causal or between-trait contrast claim.",side=1,line=9.4,adj=0,cex=.79)
  mtext("H and x are withheld interpretations, not zero effects. Original full 396-pair testing family remains unchanged.",side=1,line=10.4,adj=0,cex=.79)
}
svg(paste0(args[2],".svg"),width=15,height=9.2);draw();dev.off()
png(paste0(args[2],".png"),width=3000,height=1840,res=200,type="cairo");draw();dev.off()

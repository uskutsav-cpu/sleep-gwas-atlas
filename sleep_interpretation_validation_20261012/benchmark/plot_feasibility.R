# Plot actual aggregate counts only; no intervals or invented observations.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('Usage: Rscript plot_feasibility.R TABLE_TSV OUTPUT_PREFIX')
dat <- read.delim(args[1], stringsAsFactors = FALSE, check.names = FALSE)
models <- c('gpt-6.1-sol', 'gpt-6-sol', 'gpt-6-luna')
if (nrow(dat) != 9 || any(dat$planned != 10)) stop('Expected registered feasibility table')
draw <- function() {
  par(mfrow=c(1,2), mar=c(4.3,8.1,4.2,1.1), oma=c(3.5,0,1.5,0),
      family='sans', cex=0.92, las=1)
  b <- dat[match(paste(models, 'B'), paste(dat$model_alias, dat$condition)),]
  c <- dat[match(paste(models, 'C'), paste(dat$model_alias, dat$condition)),]
  vals <- rbind(b$strict_correct / b$completed, c$strict_correct / c$completed)
  mids <- barplot(vals, beside=TRUE, horiz=TRUE, xlim=c(0,1.25),
                  names.arg=models, col=c('#718096','#147D92'), border=NA,
                  axes=FALSE, xlab='Strict correct task fraction',
                  main='(a)   Matched evidence: B versus C')
  axis(1, at=seq(0,1,0.25), labels=c('0','0.25','0.50','0.75','1.00'))
  labs <- rbind(paste0(b$strict_correct,'/',b$completed),
                paste0(c$strict_correct,'/',c$completed))
  text(as.vector(vals)+0.055, as.vector(mids), as.vector(labs), pos=4, cex=0.86)
  legend(x=0, y=par('usr')[4]+0.55, xpd=NA, horiz=TRUE,
         legend=c('B: disjoint tables','C: joined atlas'),
         fill=c('#718096','#147D92'), bty='n', cex=0.82)
  a <- dat[match(paste(models, 'A'), paste(dat$model_alias, dat$condition)),]
  av <- a$strict_correct/a$completed
  am <- barplot(av, horiz=TRUE, xlim=c(0,1.25), names.arg=models,
                col='#9A7541', border=NA, axes=FALSE,
                xlab='Strict correct task fraction',
                main='(b)   Conventional resource: A')
  axis(1, at=seq(0,1,0.25), labels=c('0','0.25','0.50','0.75','1.00'))
  text(av+0.055, am, paste0(a$strict_correct,'/',a$completed), pos=4, cex=0.86)
  mtext('A has no authenticated exact target-source coverage; correct abstention is credited.',
        side=1, outer=TRUE, line=0.8, cex=0.84)
  mtext('90 completed cells; 88 scored (2 A provenance exclusions). Three OpenAI aliases; descriptive key agreement only.',
        side=1, outer=TRUE, line=2.1, cex=0.84)
}
svg(paste0(args[2], '.svg'), width=13.5, height=5.5, bg='white'); draw(); dev.off()
png(paste0(args[2], '.png'), width=2700, height=1100, res=200, bg='white'); draw(); dev.off()

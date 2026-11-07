#!/usr/bin/env Rscript
root<-'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay/pc_threshold/component_bands'
.libPaths(c('/Volumes/Brain6Work/lib','/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v1/Rlib',.libPaths()))
suppressPackageStartupMessages(library(jsonlite))
x<-fromJSON(file.path(root,'projections.json'));s<-x$squared_projections;sig<-x$sigma[4];stopifnot(length(s)==306,x$threshold_K[4]==216,sig>0)
band<-sum(s[143:216]);total<-sum(s[1:216]);stopifnot(abs(total-x$cumulative[4])<1e-12)
out<-data.frame(region=c('all_216','band_143_216'),k=c(216L,74L),sum_squares=c(total,band),sigma=sig)
out$chi_square<-out$sum_squares/out$sigma;out$lower_tail<-pchisq(out$chi_square,out$k);out$ratio<-out$chi_square/out$k
write.table(out,file.path(root,'idealized_null','results.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('sampling-term fractional reduction to zero:',1-out$ratio[1],'\n');print(sessionInfo())

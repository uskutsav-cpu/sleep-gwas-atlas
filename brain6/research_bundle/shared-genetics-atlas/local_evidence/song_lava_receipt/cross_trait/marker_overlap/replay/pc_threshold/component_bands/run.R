#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(LAVA)); suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
parent <- 'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay'
out <- file.path(parent,'pc_threshold/component_bands')
prior <- fromJSON(file.path(parent,'insomnia.json'))
previous <- fromJSON(file.path(parent,'pc_threshold/results.json'))
base <- '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
input <- process.input(input.info.file=file.path(parent,'insomnia_quoted_input.tsv'),
 sample.overlap.file=file.path(base,'pairs/insomnia/sample_overlap.txt'),
 ref.prefix='/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3',
 phenos=c('frailty','insomnia'),input.dir=NULL)
set.seed(prior$seed)
loc <- process.locus(data.frame(LOC='464',CHR=3L,START=47588462L,STOP=50387742L),input,
 min.K=2L,prune.thresh=99.5,max.prop.K=.75,drop.failed=FALSE,max.block.size=3000L,cap.estimates=TRUE)
stopifnot(loc$K==306,loc$n.snps==2741,setequal(loc$snps,prior$snps))
values <- unname(loc$delta[,'frailty']^2);cum <- cumsum(values)
stopifnot(all(abs(cum[previous$K]-previous$fi_delta_sumsq)<1e-12))
write_json(list(squared_projections=values,threshold_K=previous$K,cumulative=cum[previous$K],
               sigma=previous$fi_sigma,n_snps=loc$n.snps),file.path(out,'projections.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
cat('COMPLETE; prefix sums match all five previous runs\n')

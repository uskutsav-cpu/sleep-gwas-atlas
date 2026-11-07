#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(LAVA))
suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
out <- 'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay/restricted_support'
parent <- dirname(out)
apnea <- fromJSON(file.path(parent,'sleep_apnea.json'))
prior <- fromJSON(file.path(parent,'insomnia.json'))
base <- '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
overlap <- file.path(base,'pairs/insomnia/sample_overlap.txt')
prefix <- '/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3'
info <- file.path(parent,'insomnia_quoted_input.tsv')
input <- process.input(input.info.file=info,sample.overlap.file=overlap,ref.prefix=prefix,
                       phenos=c('frailty','insomnia'),input.dir=NULL)
target <- apnea$snps
stopifnot(length(target)==2452,length(unique(target))==2452,all(target %in% input$analysis.snps))
input$analysis.snps <- intersect(input$analysis.snps,target)
set.seed(prior$seed)
loc <- process.locus(data.frame(LOC='464',CHR=3L,START=47588462L,STOP=50387742L),input,
                     min.K=2L,prune.thresh=99,max.prop.K=0.75,drop.failed=FALSE,
                     max.block.size=3000L,cap.estimates=TRUE)
stopifnot(!is.null(loc),setequal(loc$snps,target),loc$n.snps==2452,
          identical(loc$phenos,c('frailty','insomnia')))
result <- list(n_snps=loc$n.snps,K=loc$K,phenotypes=loc$phenos,N=loc$N,failed=loc$failed,
 sigma=loc$sigma,omega=loc$omega,delta_squared_sums=colSums(loc$delta^2),
 h2_obs_capped=loc$h2.obs,nref_scale=loc$nref.scale,snps=loc$snps,
 note='post-selected support restriction, diagnostic only; no bivariate test')
write_json(result,file.path(out,'restricted.json'),pretty=TRUE,auto_unbox=TRUE,digits=15,na='null')
cat('RESTRICTED_DIAGNOSTIC_COMPLETED\n')

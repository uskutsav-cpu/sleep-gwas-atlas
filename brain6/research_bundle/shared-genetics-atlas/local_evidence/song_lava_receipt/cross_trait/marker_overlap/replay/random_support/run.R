#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(LAVA))
suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
out <- 'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay/random_support'
parent <- dirname(out)
subsets <- fromJSON(file.path(out,'subsets.json'),simplifyVector=FALSE)
prior <- fromJSON(file.path(parent,'insomnia.json'))
base <- '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
prefix <- '/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3'
input <- process.input(input.info.file=file.path(parent,'insomnia_quoted_input.tsv'),
 sample.overlap.file=file.path(base,'pairs/insomnia/sample_overlap.txt'),ref.prefix=prefix,
 phenos=c('frailty','insomnia'),input.dir=NULL)
initial <- input$analysis.snps
stopifnot(length(subsets)==20)
results <- list()
for(i in seq_along(subsets)) {
 target <- unlist(subsets[[i]],use.names=FALSE)
 stopifnot(length(target)==2452,length(unique(target))==2452,all(target %in% initial))
 input$analysis.snps <- intersect(initial,target)
 set.seed(prior$seed)
 loc <- process.locus(data.frame(LOC='464',CHR=3L,START=47588462L,STOP=50387742L),input,
  min.K=2L,prune.thresh=99,max.prop.K=0.75,drop.failed=FALSE,max.block.size=3000L,cap.estimates=TRUE)
 stopifnot(!is.null(loc),loc$n.snps==2452,setequal(loc$snps,target),identical(loc$phenos,c('frailty','insomnia')))
 results[[i]] <- list(replicate=i,K=loc$K,fi_omega=loc$omega['frailty','frailty'],
  fi_failed=loc$failed['frailty'],fi_h2_capped=loc$h2.obs['frailty'],
  fi_delta_squared=sum(loc$delta[,'frailty']^2))
 cat('replicate',i,'complete\n');flush.console()
}
write_json(results,file.path(out,'results.json'),pretty=TRUE,auto_unbox=TRUE,digits=15,na='null')

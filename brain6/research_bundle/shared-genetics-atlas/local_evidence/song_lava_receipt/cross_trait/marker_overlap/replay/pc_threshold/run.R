#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(LAVA))
suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
parent <- 'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay'
out <- file.path(parent,'pc_threshold')
prior <- fromJSON(file.path(parent,'insomnia.json'))
base <- '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
input <- process.input(input.info.file=file.path(parent,'insomnia_quoted_input.tsv'),
 sample.overlap.file=file.path(base,'pairs/insomnia/sample_overlap.txt'),
 ref.prefix='/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3',
 phenos=c('frailty','insomnia'),input.dir=NULL)
rows <- list()
for (threshold in c(95,97,98,99,99.5)) {
 set.seed(prior$seed)
 loc <- process.locus(data.frame(LOC='464',CHR=3L,START=47588462L,STOP=50387742L),input,
   min.K=2L,prune.thresh=threshold,max.prop.K=.75,drop.failed=FALSE,max.block.size=3000L,cap.estimates=TRUE)
 stopifnot(!is.null(loc),loc$n.snps==2741,setequal(loc$snps,prior$snps),identical(loc$phenos,c('frailty','insomnia')))
 rows[[length(rows)+1]] <- list(threshold=threshold,n_snps=loc$n.snps,K=loc$K,
  fi_omega=unname(loc$omega['frailty','frailty']),fi_sigma=unname(loc$sigma['frailty','frailty']),
  fi_h2_capped=unname(loc$h2.obs['frailty']),fi_delta_sumsq=unname(colSums(loc$delta^2)['frailty']),
  fi_failed=unname(loc$failed['frailty']))
 if(threshold==99)stopifnot(loc$K==prior$K,abs(loc$omega['frailty','frailty']-prior$omega[1,1])<1e-12)
 cat('completed',threshold,'K',loc$K,'fi_omega',loc$omega['frailty','frailty'],'\n');flush.console()
}
write_json(rows,file.path(out,'results.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
cat('COMPLETE\n')

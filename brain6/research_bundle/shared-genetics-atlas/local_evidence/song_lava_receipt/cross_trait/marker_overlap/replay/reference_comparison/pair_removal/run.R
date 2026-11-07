suppressPackageStartupMessages(library(LAVA));suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
root <- 'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay'
out <- file.path(root,'reference_comparison/pair_removal');prior <- fromJSON(file.path(root,'insomnia.json'))
base <- '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
input <- process.input(input.info.file=file.path(root,'insomnia_quoted_input.tsv'),
 sample.overlap.file=file.path(base,'pairs/insomnia/sample_overlap.txt'),
 ref.prefix='/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3',
 phenos=c('frailty','insomnia'),input.dir=NULL)
original <- input$analysis.snps
remove <- c('rs62262671','rs55754265');stopifnot(all(remove %in% prior$snps),length(prior$snps)==2741)
results <- list()
for (id in remove){
 input$analysis.snps <- setdiff(original,id);set.seed(prior$seed)
 loc<-process.locus(data.frame(LOC='464',CHR=3L,START=47588462L,STOP=50387742L),input,
  min.K=2L,prune.thresh=99,max.prop.K=.75,drop.failed=FALSE,max.block.size=3000L,cap.estimates=TRUE)
 stopifnot(!is.null(loc),loc$n.snps==2740,setequal(loc$snps,setdiff(prior$snps,id)))
 results[[length(results)+1]]<-list(removed=id,n_snps=loc$n.snps,K=loc$K,
  fi_omega=unname(loc$omega['frailty','frailty']),fi_sigma=unname(loc$sigma['frailty','frailty']),
  fi_h2_capped=unname(loc$h2.obs['frailty']),fi_delta_sumsq=unname(colSums(loc$delta^2)['frailty']),
  fi_failed=unname(loc$failed['frailty']))
 cat('completed',id,'K',loc$K,'omega',loc$omega['frailty','frailty'],'\n');flush.console()
}
write_json(results,file.path(out,'results.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
cat('PAIR_REMOVAL_COMPLETE\n')

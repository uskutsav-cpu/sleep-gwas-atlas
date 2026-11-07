suppressPackageStartupMessages(library(LAVA));suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
root <- 'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay'
out<-file.path(root,'reference_comparison/eur_lava_sensitivity');prior<-fromJSON(file.path(root,'insomnia.json'))
support<-read.delim(file.path(out,'support.tsv'),stringsAsFactors=FALSE);eur<-fromJSON(file.path(out,'eur_result.json'))
stopifnot(nrow(support)==eur$matched_markers,all(support$SNP %in% prior$snps))
base<-'/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
input<-process.input(input.info.file=file.path(root,'insomnia_quoted_input.tsv'),
 sample.overlap.file=file.path(base,'pairs/insomnia/sample_overlap.txt'),
 ref.prefix='/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3',
 phenos=c('frailty','insomnia'),input.dir=NULL)
input$analysis.snps<-intersect(input$analysis.snps,support$SNP);set.seed(prior$seed)
loc<-process.locus(data.frame(LOC='464',CHR=3L,START=47588462L,STOP=50387742L),input,
 min.K=2L,prune.thresh=99,max.prop.K=.75,drop.failed=FALSE,max.block.size=3000L,cap.estimates=TRUE)
stopifnot(!is.null(loc),identical(loc$snps,support$SNP),loc$n.snps==nrow(support))
output<-list(n_snps=loc$n.snps,K=loc$K,N_fi=unname(loc$N['frailty']),
 fi_omega=unname(loc$omega['frailty','frailty']),fi_sigma=unname(loc$sigma['frailty','frailty']),
 fi_h2_capped=unname(loc$h2.obs['frailty']),fi_delta_sumsq=unname(colSums(loc$delta^2)['frailty']),fi_failed=unname(loc$failed['frailty']))
write_json(output,file.path(out,'matched_ukb.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
cat('MATCHED_UKB_COMPLETE\n')

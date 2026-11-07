#!/usr/bin/env Rscript
# Diagnostic replay only: retain failed phenotypes; never run bivariate tests.
suppressPackageStartupMessages(library(LAVA))
suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion()) == '4.3.3', as.character(packageVersion('LAVA')) == '0.1.5')
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1, args[[1]] %in% c('insomnia','sleep_apnea'))
trait <- args[[1]]
base <- '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
pair <- file.path(base,'pairs',trait)
original_info <- file.path(pair,'input_info_chr03.tsv')
info <- file.path('outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay',paste0(trait,'_quoted_input.tsv'))
# read.table in LAVA splits unquoted paths with spaces. Preserve the saved fields
# but quote all character fields and use the saved absolute paths directly.
descriptor <- read.delim(original_info,stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(identical(descriptor$phenotype,c('frailty',trait)),all(file.exists(descriptor$filename)))
write.table(descriptor,info,sep='\t',row.names=FALSE,quote=TRUE,na='NA')
overlap <- file.path(pair,'sample_overlap.txt')
prefix <- '/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3'
locrow <- data.frame(LOC='464',CHR=3L,START=47588462L,STOP=50387742L)
saved <- fromJSON(file.path(pair,'loci/locus_0464.json'))
set.seed(saved$seed)
input <- process.input(input.info.file=info,sample.overlap.file=overlap,
                       ref.prefix=prefix,phenos=c('frailty',trait),input.dir=NULL)
loc <- process.locus(locrow,input,min.K=2L,prune.thresh=99,max.prop.K=0.75,
                     drop.failed=FALSE,max.block.size=3000L,cap.estimates=TRUE)
stopifnot(!is.null(loc),loc$n.snps==saved$n_snps,loc$K==saved$K,
          identical(loc$phenos,c('frailty',trait)))
result <- list(trait=trait, seed=saved$seed, n_snps=loc$n.snps, K=loc$K,
  phenotypes=loc$phenos,N=loc$N,failed=loc$failed,h2_obs_capped=loc$h2.obs,
  sigma=loc$sigma,omega=loc$omega,nref_scale=loc$nref.scale,
  delta_squared_sums=colSums(loc$delta^2),snps=loc$snps,
  historical_sleep_h2=saved$sleep_h2_obs,historical_fi_h2=saved$fi_h2_obs,
  historical_sleep_status=saved$sleep_status,historical_fi_status=saved$fi_status,
  warning='Diagnostic replay, not a bivariate test or a gate rescue')
stopifnot(abs(loc$h2.obs[trait]-saved$sleep_h2_obs) < 1e-9)
if (trait=='sleep_apnea') stopifnot(abs(loc$h2.obs['frailty']-saved$fi_h2_obs)<1e-9)
out <- file.path('outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay',paste0(trait,'.json'))
write_json(result,out,pretty=TRUE,auto_unbox=TRUE,digits=15,na='null')
cat('REPLAY_COMPLETED',trait,'\n')

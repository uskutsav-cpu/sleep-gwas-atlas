suppressPackageStartupMessages(library(LAVA));suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
root <- 'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay'
out <- file.path(root,'reference_comparison'); prior <- fromJSON(file.path(root,'insomnia.json'))
base <- '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
input <- process.input(input.info.file=file.path(root,'insomnia_quoted_input.tsv'),
 sample.overlap.file=file.path(base,'pairs/insomnia/sample_overlap.txt'),
 ref.prefix='/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3',
 phenos=c('frailty','insomnia'),input.dir=NULL)
idx <- round(seq(1,length(prior$snps),length.out=60));stopifnot(length(unique(idx))==60)
ld <- LAVA:::read.ld(input$reference,prior$snps[idx]);stopifnot(nrow(ld$info)==60,identical(ld$info$SNP,prior$snps[idx]))
write.table(ld$info[,c('SNP','CHR','POS','A1','A2','FREQ')],file.path(out,'markers.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(SNP=rownames(ld$ld),ld$ld,check.names=FALSE),file.path(out,'ukb_ld.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('UKB_LD_COMPLETE\n')

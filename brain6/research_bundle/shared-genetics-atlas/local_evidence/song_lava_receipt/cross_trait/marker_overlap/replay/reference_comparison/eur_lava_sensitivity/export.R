suppressPackageStartupMessages(library(LAVA));suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
root <- 'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay'
out<-file.path(root,'reference_comparison/eur_lava_sensitivity');prior<-fromJSON(file.path(root,'insomnia.json'))
base<-'/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
input<-process.input(input.info.file=file.path(root,'insomnia_quoted_input.tsv'),
 sample.overlap.file=file.path(base,'pairs/insomnia/sample_overlap.txt'),
 ref.prefix='/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3',
 phenos=c('frailty','insomnia'),input.dir=NULL)
ref<-input$reference$snp.info;ref<-ref[match(prior$snps,ref$SNP),]
fi<-input$sum.stats[['frailty']];fi<-fi[match(prior$snps,fi$SNP),]
stopifnot(identical(ref$SNP,prior$snps),identical(fi$SNP,prior$snps),all(c('STAT','N') %in% names(fi)))
raw<-LAVA:::read.sumstats.file(input$info$filename[input$info$phenotype=='frailty'],'frailty')
raw<-raw[match(prior$snps,raw$SNP),]
orient<-LAVA:::com.pair(LAVA:::map.alleles(ref$A1,ref$A2),LAVA:::map.alleles(raw$A1,raw$A2))
stopifnot(all(is.finite(orient)),all(abs(fi$STAT - raw$STAT*orient)<1e-12))
result<-data.frame(SNP=ref$SNP,CHR=ref$CHR,POS=ref$POS,A1=ref$A1,A2=ref$A2,STAT=fi$STAT,N=fi$N)
stopifnot(all(is.finite(result$STAT)),all(is.finite(result$N)))
write.table(result,file.path(out,'aligned_fi.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('ALIGNED_FI_EXPORTED\n')

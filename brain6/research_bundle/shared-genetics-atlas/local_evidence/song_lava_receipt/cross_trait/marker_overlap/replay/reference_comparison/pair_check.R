suppressPackageStartupMessages(library(LAVA));suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
root <- 'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay'
base <- '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
input <- process.input(input.info.file=file.path(root,'insomnia_quoted_input.tsv'),
 sample.overlap.file=file.path(base,'pairs/insomnia/sample_overlap.txt'),
 ref.prefix='/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3',
 phenos=c('frailty','insomnia'),input.dir=NULL)
ids<-c('rs62262671','rs55754265');m<-LAVA:::read.ld(input$reference,ids)$ld
stopifnot(identical(rownames(m),ids),abs(m[1,2]-1.0016932024945)<1e-12)
write_json(list(ids=ids,matrix=unname(m),offdiagonal=unname(m[1,2])),file.path(root,'reference_comparison/pair_check.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
cat('PAIR_REPRODUCED\n')

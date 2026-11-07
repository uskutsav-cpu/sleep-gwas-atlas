suppressPackageStartupMessages(library(LAVA));suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
root <- 'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay'
out <- file.path(root,'reference_comparison/full_spectrum');prior<-fromJSON(file.path(root,'insomnia.json'))
base <- '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
input <- process.input(input.info.file=file.path(root,'insomnia_quoted_input.tsv'),
 sample.overlap.file=file.path(base,'pairs/insomnia/sample_overlap.txt'),
 ref.prefix='/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3',
 phenos=c('frailty','insomnia'),input.dir=NULL)
ld<-LAVA:::read.ld(input$reference,prior$snps);m<-ld$ld
stopifnot(nrow(m)==2741,identical(rownames(m),prior$snps),all(is.finite(m)))
asym<-max(abs(m-t(m)));original_diag<-diag(m);diagdev<-max(abs(original_diag-1));diag(m)<-0
aboveone<-sum(abs(m)>1)/2;maxabs<-max(abs(m));diag(m)<-original_diag
cat('matrix checked; computing spectrum\n');flush.console()
ev<-eigen(m,symmetric=TRUE,only.values=TRUE)$values
cum<-cumsum(ev)/max(cumsum(ev))*100;k<-min(which(cum>=99))
stopifnot(k==prior$K,all(ev[1:k]>0))
result<-list(n_snps=nrow(m),max_asymmetry=asym,max_diagonal_deviation=diagdev,
 offdiagonal_abs_gt_1=aboveone,max_abs_offdiagonal=maxabs,
 minimum_eigenvalue=min(ev),negative_below_1e_8=sum(ev < -1e-8),negative_below_1e_4=sum(ev < -1e-4),
 threshold_99_K=k,smallest_retained_eigenvalue=ev[k],eigenvalues=ev,
 note='Approximate reference matrix; no local bivariate correlation or quality-gate repair.')
write_json(result,file.path(out,'spectrum.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
cat('FULL_SPECTRUM_COMPLETE\n')

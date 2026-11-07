suppressPackageStartupMessages(library(dplyr));suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(jsonlite))
work<-'/Volumes/Brain6Work/hdll-chr3-44-20261028';out<-'outputs/shared-genetics-atlas/local_evidence/hdll_pilot'
source(file.path(work,'HDL.L.R'))
base<-'/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1/inputs'
ref<-read.table(file.path(work,'bim/ukb_chr3.44_n336000.imputed_clean.bim'),stringsAsFactors=FALSE)
colnames(ref)<-c('CHR','SNP','CM','POS','A1','A2');stopifnot(nrow(ref)==681,!anyDuplicated(ref$SNP))
readtrait<-function(t){x<-as.data.frame(fread(file.path(base,t,'chr03.sumstats.tsv.gz')));x<-x[x$SNP %in% ref$SNP,c('SNP','A1','A2','Z','N')];stopifnot(!anyDuplicated(x$SNP),all(is.finite(x$Z)),all(is.finite(x$N)));x}
fi<-readtrait('frailty');ins<-readtrait('insomnia')
check<-function(x){z<-ref[match(x$SNP,ref$SNP),];all((x$A1==z$A1 & x$A2==z$A2)|(x$A1==z$A2 & x$A2==z$A1))}
stopifnot(nrow(fi)==678,nrow(ins)==677,check(fi),check(ins))
results<-list();logs<-list()
for(n0 in c(0,50000,100000)){
 warnings<-character();value<-NULL;err<-NULL
 logs[[as.character(n0)]]<-capture.output(value<-tryCatch(withCallingHandlers(
  HDL.L(fi,ins,Trait1name='FI',Trait2name='insomnia',LD.path=file.path(work,'ld'),bim.path=file.path(work,'bim'),
   Nref=335272,N0=n0,chr=3,piece=44,output.file='',eigen.cut=.99),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
  error=function(e){err<<-conditionMessage(e);NULL}))
 results[[length(results)+1]]<-list(N0_assumed=n0,status=if(is.null(err))'returned' else 'error',error=err,warnings=warnings,value=value)
 cat('completed N0',n0,'status',if(is.null(err))'returned' else err,'\n');flush.console()
}
write_json(list(reference_snps=681,FI_matched=nrow(fi),insomnia_matched=nrow(ins),
 FI_median_N=median(fi$N),insomnia_median_N=median(ins$N),runs=results,logs=logs,
 note='Post-selected sensitivity, hypothetical overlap; no validated inference or QC gate rescue.'),
 file.path(out,'run_results.json'),pretty=TRUE,auto_unbox=TRUE,digits=15,na='null')
cat('HDLL_PILOT_COMPLETE\n')

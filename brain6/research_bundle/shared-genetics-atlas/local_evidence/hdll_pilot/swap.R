suppressPackageStartupMessages(library(dplyr));suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(jsonlite))
work<-'/Volumes/Brain6Work/hdll-chr3-44-20261028';out<-'outputs/shared-genetics-atlas/local_evidence/hdll_pilot'
source(file.path(work,'HDL.L.R'));old<-fromJSON(file.path(out,'run_results.json'),simplifyVector=FALSE)
ref<-read.table(file.path(work,'bim/ukb_chr3.44_n336000.imputed_clean.bim'),stringsAsFactors=FALSE);ids<-ref[,2]
base<-'/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1/inputs'
readtrait<-function(t){x<-as.data.frame(fread(file.path(base,t,'chr03.sumstats.tsv.gz')));x[x$SNP %in% ids,c('SNP','A1','A2','Z','N')]}
fi<-readtrait('frailty');ins<-readtrait('insomnia');warnings<-character();v<-NULL;err<-NULL
log<-capture.output(v<-tryCatch(withCallingHandlers(HDL.L(ins,fi,Trait1name='insomnia',Trait2name='FI',LD.path=file.path(work,'ld'),bim.path=file.path(work,'bim'),Nref=335272,N0=50000,chr=3,piece=44,output.file='',eigen.cut=.99),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);NULL}))
write_json(list(status=if(is.null(err))'returned' else 'error',error=err,warnings=warnings,value=v,log=log),file.path(out,'swap_results.json'),pretty=TRUE,auto_unbox=TRUE,digits=15,na='null')
cat('SWAP_CHECK_COMPLETE\n')

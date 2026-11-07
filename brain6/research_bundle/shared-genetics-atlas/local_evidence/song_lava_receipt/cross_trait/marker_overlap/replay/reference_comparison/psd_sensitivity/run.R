suppressPackageStartupMessages(library(LAVA));suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
root <- 'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay'
out <- file.path(root,'reference_comparison/psd_sensitivity');prior<-fromJSON(file.path(root,'insomnia.json'))
base <- '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
input <- process.input(input.info.file=file.path(root,'insomnia_quoted_input.tsv'),
 sample.overlap.file=file.path(base,'pairs/insomnia/sample_overlap.txt'),
 ref.prefix='/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3',
 phenos=c('frailty','insomnia'),input.dir=NULL)
raw <- LAVA:::read.ld(input$reference,prior$snps);m<-raw$ld
stopifnot(nrow(m)==2741,identical(rownames(m),prior$snps),max(abs(m-t(m)))<1e-12)
orig_read <- getFromNamespace('read.ld','LAVA')
wrapper <- function(reference,snp.ids,require.freq=FALSE){
 stopifnot(setequal(snp.ids,prior$snps)); list(mode='ld',info=raw$info,ld=current_matrix)
}
assignInNamespace('read.ld',wrapper,ns='LAVA')
normalise <- function(x){s<-sqrt(diag(x));stopifnot(all(s>0));x/outer(s,s)}
results<-list();ridge<-NA_real_
for(setting in c('unchanged','clipped','ridge')){
 if(setting=='unchanged')current_matrix<-m
 if(setting=='clipped'){
  e<-eigen(m,symmetric=TRUE);ridge<-max(0,-min(e$values))+1e-8
  current_matrix<-normalise((e$vectors * rep(pmax(e$values,0),each=nrow(m))) %*% t(e$vectors))
  rm(e);gc()
 }
 if(setting=='ridge'){current_matrix<-m;diag(current_matrix)<-diag(current_matrix)+ridge;current_matrix<-normalise(current_matrix)}
 dimnames(current_matrix)<-dimnames(m)
 change<-current_matrix-m;max_change<-max(abs(change));rel_frob<-sqrt(sum(change^2)/sum(m^2));rm(change);gc()
 minimum<-if(setting=='unchanged')NA_real_ else min(eigen(current_matrix,symmetric=TRUE,only.values=TRUE)$values)
 if(setting!='unchanged')stopifnot(minimum > -1e-7,max(abs(diag(current_matrix)-1))<1e-12)
 set.seed(prior$seed)
 loc<-process.locus(data.frame(LOC='464',CHR=3L,START=47588462L,STOP=50387742L),input,
  min.K=2L,prune.thresh=99,max.prop.K=.75,drop.failed=FALSE,max.block.size=3000L,cap.estimates=TRUE)
 stopifnot(loc$n.snps==2741,identical(loc$snps,prior$snps))
 if(setting=='unchanged')stopifnot(loc$K==prior$K,abs(loc$omega['frailty','frailty']-prior$omega[1,1])<1e-12)
 results[[length(results)+1]]<-list(setting=setting,minimum_eigenvalue=minimum,ridge=if(setting=='ridge')ridge else NA_real_,
  relative_frobenius_change=rel_frob,max_entry_change=max_change,n_snps=loc$n.snps,K=loc$K,
  fi_omega=unname(loc$omega['frailty','frailty']),fi_sigma=unname(loc$sigma['frailty','frailty']),
  fi_h2_capped=unname(loc$h2.obs['frailty']),fi_delta_sumsq=unname(colSums(loc$delta^2)['frailty']),fi_failed=unname(loc$failed['frailty']))
 cat('completed',setting,'K',loc$K,'omega',loc$omega['frailty','frailty'],'\n');flush.console()
}
assignInNamespace('read.ld',orig_read,ns='LAVA')
write_json(results,file.path(out,'results.json'),pretty=TRUE,auto_unbox=TRUE,digits=15,na='null')
cat('PSD_SENSITIVITY_COMPLETE\n')

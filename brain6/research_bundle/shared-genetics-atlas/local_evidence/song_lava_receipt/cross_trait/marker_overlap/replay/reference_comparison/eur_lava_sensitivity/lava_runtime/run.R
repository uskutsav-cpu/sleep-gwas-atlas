suppressPackageStartupMessages(library(LAVA));suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('LAVA'))=='0.1.5')
root<-'outputs/shared-genetics-atlas/local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay'
out<-file.path(root,'reference_comparison/eur_lava_sensitivity/lava_runtime');parent<-dirname(out)
support<-read.delim(file.path(parent,'support.tsv'),stringsAsFactors=FALSE);prior<-fromJSON(file.path(root,'insomnia.json'));eur<-fromJSON(file.path(parent,'eur_result.json'))
n<-nrow(support);stopifnot(n==2738)
con<-file('stdin','rb');raw<-readBin(con,what='double',n=n*n,size=8,endian='little');close(con)
stopifnot(length(raw)==n*n);ldmat<-matrix(raw,nrow=n,ncol=n);rm(raw)
stopifnot(max(abs(ldmat-t(ldmat)))<1e-12,max(abs(diag(ldmat)-1))<1e-10)
dimnames(ldmat)<-list(support$SNP,support$SNP)
base<-'/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1'
input<-process.input(input.info.file=file.path(root,'insomnia_quoted_input.tsv'),
 sample.overlap.file=file.path(base,'pairs/insomnia/sample_overlap.txt'),
 ref.prefix='/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3',
 phenos=c('frailty','insomnia'),input.dir=NULL)
input$analysis.snps<-intersect(input$analysis.snps,support$SNP)
info<-input$reference$snp.info[match(support$SNP,input$reference$snp.info$SNP),]
stopifnot(identical(info$SNP,support$SNP))
orig<-getFromNamespace('read.ld','LAVA')
wrapper<-function(reference,snp.ids,require.freq=FALSE){stopifnot(identical(snp.ids,support$SNP));list(mode='ld',info=info,ld=ldmat)}
assignInNamespace('read.ld',wrapper,ns='LAVA');set.seed(prior$seed)
loc<-process.locus(data.frame(LOC='464',CHR=3L,START=47588462L,STOP=50387742L),input,
 min.K=2L,prune.thresh=99,max.prop.K=.75,drop.failed=FALSE,max.block.size=3000L,cap.estimates=TRUE)
assignInNamespace('read.ld',orig,ns='LAVA')
stopifnot(loc$n.snps==n,loc$K==eur$K_99,abs(loc$omega['frailty','frailty']-eur$fi_omega)<1e-12)
ok<-loc$phenos[!loc$failed];univ<-if(length(ok))run.univ(loc,phenos=ok,var=TRUE) else NULL
result<-list(n_snps=n,K=loc$K,phenotypes=loc$phenos,failed=loc$failed,N=loc$N,
 omega=loc$omega,sigma=loc$sigma,h2_obs=loc$h2.obs,delta_squared_sums=colSums(loc$delta^2),univ=univ,
 note='1000 Genomes dosage LD, original UKB reference frequencies; univariate diagnostic only')
write_json(result,file.path(out,'lava_result.json'),pretty=TRUE,auto_unbox=TRUE,digits=15,na='null')
cat('LAVA_UNIVARIATE_CHECK_COMPLETE\n')

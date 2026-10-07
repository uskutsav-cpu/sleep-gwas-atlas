lib<-normalizePath('work/R-scientific-lib',mustWork=FALSE)
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
paths<-readLines('brain6/translational_psychiatry_research_v1/research_v1/ld/R_install_source_order.txt')
for(path in paths) {
 cat('INSTALLING_EXACT_SOURCE ',basename(path),'\n')
 install.packages(path,repos=NULL,type='source',lib=lib,Ncpus=1)
 size<-sum(file.info(list.files(lib,recursive=TRUE,full.names=TRUE))$size,na.rm=TRUE)
 if(size>150*1024^2)stop('Installed library exceeded150MiBcap')
}
stopifnot(as.character(packageVersion('susieR'))=='0.14.2',as.character(packageVersion('coloc'))=='5.2.3')
writeLines(capture.output(sessionInfo()),'brain6/translational_psychiatry_research_v1/research_v1/ld/R_runtime_sessionInfo.txt')
write.table(installed.packages(lib.loc=lib)[,c('Package','Version','Built'),drop=FALSE],'brain6/translational_psychiatry_research_v1/research_v1/ld/R_installed_packages.tsv',sep='\t',quote=FALSE,row.names=FALSE)

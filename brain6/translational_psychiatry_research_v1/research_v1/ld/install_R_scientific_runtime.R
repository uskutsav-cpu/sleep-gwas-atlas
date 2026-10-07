args <- commandArgs(trailingOnly=TRUE)
lib <- normalizePath(args[1],mustWork=FALSE)
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
avail <- read.dcf(gzfile('work/ld_genotypes_research_v1/CRAN_SOURCE_PACKAGES.gz'))
rownames(avail) <- avail[,'Package']
need <- unique(c('susieR','coloc',unlist(tools::package_dependencies(c('susieR','coloc'),db=avail,which=c('Depends','Imports','LinkingTo'),recursive=TRUE))))
installed <- rownames(installed.packages())
need <- setdiff(need,c(installed,'R'))
fields <- intersect(c('Package','Version','MD5sum','Depends','Imports','LinkingTo','NeedsCompilation'),colnames(avail))
clean <- avail[need,fields,drop=FALSE]
clean[] <- gsub('\\n',' ',clean)
write.table(clean,'brain6/translational_psychiatry_research_v1/research_v1/ld/R_dependency_manifest.tsv',sep='\t',quote=FALSE,row.names=FALSE)
print(need)
if (length(args)>1 && args[2]=='install') {
 options(timeout=120,download.file.method='libcurl',Ncpus=1)
 install.packages(need,repos='https://cran.r-project.org',lib=lib,type='source',dependencies=FALSE)
 stopifnot(as.character(packageVersion('susieR'))=='0.14.2',as.character(packageVersion('coloc'))=='5.2.3')
 writeLines(capture.output(sessionInfo()),'brain6/translational_psychiatry_research_v1/research_v1/ld/R_runtime_sessionInfo.txt')
 write.table(installed.packages(lib.loc=lib)[,c('Package','Version','Built'),drop=FALSE],'brain6/translational_psychiatry_research_v1/research_v1/ld/R_installed_packages.tsv',sep='\t',quote=FALSE,row.names=FALSE)
}

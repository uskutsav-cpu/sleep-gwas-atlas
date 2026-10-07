.libPaths(c(normalizePath('work/R-scientific-lib'),.libPaths()))
local<-'work/ld_genotypes_research_v1/C_finemap_inputs';models<-c('PRIMARY_L10_MEDIAN_NEFF','SENSITIVITY_L1_MEDIAN_NEFF','SENSITIVITY_L5_MEDIAN_NEFF','SENSITIVITY_L10_N_INFINITY')
for(model in models) {
 a<-readRDS(file.path(local,paste0('insomnia_',model,'_native_fit.rds')));b<-readRDS(file.path(local,paste0('adhd_',model,'_native_fit.rds')))
 stopifnot(length(a$sets$cs)==1,length(b$sets$cs)==1)
 d<-data.frame(variant_key=colnames(a$lbf_variable),logBF1=a$lbf_variable[a$sets$cs_index[1],],logBF2=b$lbf_variable[b$sets$cs_index[1],])
 write.table(d,file.path(local,paste0(model,'_native_BFs_LOCAL_ONLY.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
}

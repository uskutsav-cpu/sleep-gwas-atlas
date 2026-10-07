.libPaths(c(normalizePath('work/R-scientific-lib'),.libPaths()))
local<-'work/ld_genotypes_research_v1/C_finemap_inputs';here<-'brain6/translational_psychiatry_research_v1/research_v1/ld'
d<-read.delim(file.path(local,'insomnia_C_shared.tsv'));p<-nrow(d)
R<-matrix(readBin(file.path(local,'shared_signed_LD.float64.bin'),what=numeric(),n=p*p,size=8,endian='little'),nrow=p)
# Deliberate reproduction of the erroneous API call; never used as a QC/pass result.
warnings<-character();error<-''
value<-tryCatch(withCallingHandlers(susieR::estimate_s_rss(z=d$Z_ALT,R=R,n=Inf,method='null-mle'),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){error<<-conditionMessage(e);NA_real_})
write.table(data.frame(warning_index=seq_along(warnings),message=warnings,deliberate_invalid_call='n=Inf',used_for_scientific_result=FALSE),file.path(here,'initial_invalid_N_infinity_reproduced_warnings.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('DELIBERATE_INVALID_INFINITY_WARNING_REPRODUCTION_COUNT ',length(warnings),' ERROR ',error,'\n')
cat('INVALID_VALUE_EXCLUDED_FROM_ALL_SCIENTIFIC_QC_OR_POSTERIOR_ANALYSES\n')

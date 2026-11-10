args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('Usage: Rscript stage5_siraj_best_library_extract.R SOURCE_RDS OUTPUT_CSV')
q <- read.csv('outputs/stage5_siraj2026_neural_mpra_exact_alleles.csv',stringsAsFactors=FALSE)
ids <- unique(paste0('chr',q$variant_GRCh37_source))
d <- readRDS(args[1])
required <- c('variant','cell_type','library','best_library','type','log2FC','logPadj_BF','Log2Skew','Skew_logPadj')
stopifnot(all(required %in% names(d)))
keep <- !is.na(d$variant) & d$variant %in% ids & !is.na(d$cell_type) & grepl('SK.?N.?SH',d$cell_type,ignore.case=TRUE) & !is.na(d$best_library) & d$best_library==1 & !is.na(d$type) & d$type!='other_test'
cols <- intersect(c(required,'variant_hg38','gene','tissue','trait','cohort','pip','ID','A_log2FC','A_logPadj_BF','B_log2FC','B_logPadj_BF','A_Ctrl_Mean','B_Ctrl_Mean','mean_Plasmid_ref','mean_Plasmid_alt','mean_RNA_ref','mean_RNA_alt'),names(d))
out <- d[keep,cols,drop=FALSE]
out$active_screen <- !is.na(out$logPadj_BF) & !is.na(out$log2FC) & out$logPadj_BF >= 2 & abs(out$log2FC) >= 1
out$emvar_screen <- out$active_screen & !is.na(out$Skew_logPadj) & out$Skew_logPadj >= 1
write.csv(out,args[2],row.names=FALSE,na='')
cat('Input rows',nrow(d),'matched rows',nrow(out),'variants',length(unique(out$variant)),'emvar rows',sum(out$emvar_screen),'emvar variants',length(unique(out$variant[out$emvar_screen])),'\n')

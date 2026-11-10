args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('Usage: Rscript stage5_siraj_all_cell_selected_extract.R SOURCE_RDS OUTPUT_CSV')
snv<-read.csv('outputs/stage5_reciprocal_snv_hg19_reference_check.csv',stringsAsFactors=FALSE)
ind<-read.csv('outputs/stage5_indel_hg19_mapping.csv',stringsAsFactors=FALSE)
ids<-unique(c(paste0('chr',snv$chr37[snv$reference_status=='reference_same'],':',snv$pos37[snv$reference_status=='reference_same'],':',snv$ref38[snv$reference_status=='reference_same'],':',snv$alt38[snv$reference_status=='reference_same']),paste0('chr',ind$hg19_chrom[ind$status=='mapped_reference_checked_normalized'],':',ind$hg19_position[ind$status=='mapped_reference_checked_normalized'],':',ind$hg19_ref[ind$status=='mapped_reference_checked_normalized'],':',ind$hg19_alt[ind$status=='mapped_reference_checked_normalized'])))
stopifnot(length(ids)==433)
d<-readRDS(args[1]);required<-c('variant','cell_type','library','best_library','type','log2FC','logPadj_BF','Log2Skew','Skew_logPadj');stopifnot(all(required %in% names(d)))
keep<-!is.na(d$variant) & d$variant %in% ids & !is.na(d$cell_type) & !is.na(d$best_library) & d$best_library==1 & !is.na(d$type) & d$type!='other_test'
cols<-intersect(c(required,'variant_hg38','gene','tissue','trait','cohort','pip','ID','A_log2FC','A_logPadj_BF','B_log2FC','B_logPadj_BF','mean_Plasmid_ref','mean_Plasmid_alt','mean_RNA_ref','mean_RNA_alt'),names(d))
out<-d[keep,cols,drop=FALSE]
out$active_screen<-!is.na(out$logPadj_BF)&!is.na(out$log2FC)&out$logPadj_BF>=2&abs(out$log2FC)>=1
out$emvar_screen<-out$active_screen & !is.na(out$Skew_logPadj)&out$Skew_logPadj>=1
write.csv(out,args[2],row.names=FALSE,na='')
cat('source rows',nrow(d),'matched rows',nrow(out),'variants',length(unique(out$variant)),'cell types',paste(unique(out$cell_type),collapse=','),'emvar variants',length(unique(out$variant[out$emvar_screen])),'\n')

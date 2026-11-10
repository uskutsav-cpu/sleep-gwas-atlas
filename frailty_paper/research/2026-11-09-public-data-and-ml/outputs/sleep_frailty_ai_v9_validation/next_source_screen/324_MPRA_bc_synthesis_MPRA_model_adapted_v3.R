
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("reshape2", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("DESeq2", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)

fileDate <- function(){
  a <- Sys.Date()
  year <- format(a, "%Y")
  month <- format(a, "%m")
  day <- format(a, "%d")
  today <- paste0(year, month, day)
  return(today)
}

### Test oligoIsolate function
### Remove Error, CIGAR, MD and position columns if necessary; aggregrate cound data with relation to the oligo
## INPUT:
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
## OUTPUT:
# oligo count data with the oligo as row names and the aggregate count data

oligoIsolate <- function(countsData){
  if("Error" %in% colnames(countsData)){
    countsData <- countsData[,c(1,2,7:dim(countsData)[2])]
  }
  tag_counts <- aggregate(. ~Oligo, data=countsData[,-1], FUN = sum)
  counts_oligo <- tag_counts[,-1]
  rownames(counts_oligo) <- tag_counts[,1]
  
  # write.table(counts_oligo, paste0(file_prefix, "_", fileDate(), "_counts.out"), quote = F, sep="\t")
  return(counts_oligo)
}

### Standardize condition data
## INPUTS:
# conditionData  : table of conditions, 2 columns no header align to the variable column headers of countsData
# and 'RNA'
# data should initially be read in such that samples are row names, if reading in the condition table produced by MPRAcount
# simply set row.names=1 when reading in the file.
## OUTPUT: standardized condition data, factorizing the cell types and categorizing DNA as 1 and RNA as 0
conditionStandard <- function(conditionData){
  cond_data <- as.data.frame(conditionData)
  colnames(cond_data)[1] <- "condition"
  cond_data[,1] <- factor(cond_data[,1])
  cond_data$condition <- relevel(cond_data$condition, "DNA")
  
  ## Not entirely sure this is necessary to keep in
  cond_data$DNA <- 0
  cond_data[cond_data$condition=="DNA",]$DNA <- 1
  ##
  
  return(cond_data)
}


### Initial processing of files via DESeq analysis
## INPUT:
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
# conditionData   : table of conditions, 2 columns no header align to the variable column headers of countsData
# and 'RNA'
# exclList        : list of celltypes that should be excluded from the analysis. Empty by default
## OUTPUT: dds_results (initial)
processAnalysis <- function(countsData, conditionData, exclList=c()){
  
  # bring in count data and condition data from previously defined functions
  
  count_data <- oligoIsolate(countsData)
  
  cat("count_data_In_function_0\n")
  cat(str(count_data))
  cat("\n")
  
  message("Oligos isolated")
  
  cond_data <- conditionStandard(conditionData)
  
  cat("cond_data_In_function_0\n")
  cat(str(cond_data))
  cat("\n")
  cat(sprintf(as.character(row.names(cond_data))))
  cat("\n")
  
  # make sure that the column names and row names of the count and condition data match
  colnames(count_data) <- row.names(cond_data)
  
  cat("count_data_In_function_1\n")
  cat(str(count_data))
  cat("\n")
  
  # perform DESeq analysis
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = cond_data, design = ~condition)
  
  message("DESeqDataSetFromMatrix")
  
  # cat("dds_In_function_0\n")
  # cat(str(dds))
  # cat("\n")
  
  dds$condition <- relevel(dds$condition, "DNA")
  dds_results <- DESeq(dds, fitType = 'local', minReplicatesForReplace=Inf)
  
  message("dds_results")
  
  
  # cat("dds_results_In_function_0\n")
  # cat(str(dds_results))
  # cat("\n")
  
  return(list(dds,dds_results))
}

### Normalize DESeq results and plot normalized densities for each cell type
## INPUT:
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
# attributesData  : table of full attributes, columns should include: ID, SNP, Project, Window, Strand, Allele,
# Haplotype, Bash
# conditionData   : table of conditions, 2 columns no header align to the variable column headers of countsData
# and 'RNA'
# exclList        : list of celltypes that should be excluded from the anlalysis. Empty by default
# method          : normalization method. Default to 'ss'
# 'ss' : Summit Shift - shifts the l2fc density peak to line up with 0
# 'ssn': Summit Shift (Negative Controls Only) - shifts the peak of negative controls to 0
# 'ro' : Remove outliers - remove oligos that don't have a p-value or have a p-value > 0.001
# 'nc' : Negative Controls - normalize only the negative controls
## OUTPUT: dds_results (normalized), plots normalization curves
tagNorm <- function(countsData, conditionData, attributesData, exclList = c(), method = 'ss', negCtrlName="negCtrl", upDisp=T, prior=F){
  process <- processAnalysis(countsData, conditionData, exclList)
  
  # cat("process_In_function_0\n")
  # cat(str(process))
  # cat("\n")
  
  
  dds_results <- process[[2]]
  
  # cat("dds_results_In_function_0\n")
  # cat(str(dds_results))
  # cat("\n")
  
  dds <- process[[1]]
  
  # cat("dds_In_function_0\n")
  # cat(str(dds))
  # cat("\n")
  
  dds_results_orig <- dds_results
  
  
  attribute_ids <- (attributesData[attributesData$ctrl_exp==negCtrlName,])$ID
  
  # cat("attribute_ids_In_function_0\n")
  # cat(str(attribute_ids))
  # cat("\n")
  
  full_attribute_ids <- attributesData$ID
  
  # cat("full_attribute_ids_In_function_0\n")
  # cat(str(full_attribute_ids))
  # cat("\n")
  
  count_data <- oligoIsolate(countsData)
  cond_data <- conditionStandard(conditionData)
  colnames(count_data) <- row.names(cond_data)
  
  # cat("dds_results_In_function_0\n")
  # cat(str(dds_results))
  # cat("\n")
  
  temp_outputA <- results(dds_results, contrast = c("condition",'RNA',"DNA"), cooksCutoff=F, independentFiltering=F)
  
  # cat("temp_outputA_In_function_0\n")
  # cat(str(temp_outputA))
  # cat("\n")
  
  # Summit shift normalization
  if(method == "ss"){
    message("summit shift")
    summit <- which.max(density(temp_outputA$log2FoldChange, na.rm=T)$y)
    log_offset <- 2^(density(temp_outputA$log2FoldChange, na.rm=T)$x[summit])
    sizeFactors(dds_results)[which(cond_data$condition == 'RNA')] <- sizeFactors(dds_results)[which(cond_data$condition == 'RNA')]*(log_offset)
  }
  # Summit shift normalization - negative controls only
  if(method == "ssn"){
    temp_outputA_neg <- temp_outputA[attribute_ids,]
    summit <- which.max(density(temp_outputA_neg$log2FoldChange, na.rm=T)$y)
    log_offset <- 2^(density(temp_outputA_neg$log2FoldChange, na.rm=T)$x[summit])
    sizeFactors(dds_results)[which(cond_data$condition == 'RNA')] <- sizeFactors(dds_results)[which(cond_data$condition == 'RNA')]*(log_offset)
  }
  # Remove outliers for normalization
  if(method == "ro"){
    temp_outputA_neg <- temp_outputA[full_attribute_ids,]
    attribute_ids <- row.names(temp_outputA_neg[!is.na(temp_outputA_neg$pvalue) & temp_outputA_neg$pvalue>0.001,])
  }
  cat("temp_outputA_In_function_1\n")
  cat(str(temp_outputA))
  cat("\n")
  
  
  
  # Remove outliers or negative controls only
  if(method == "ro" | method == "nc"){
    dds_results_tmp<-estimateSizeFactors(dds[attribute_ids])
    sizeFactors(dds_results)<-sizeFactors(dds_results_tmp)
  }
  
  message("normalization performed")
  
  return(list(dds_results,dds_results_orig))
}

tagSig <- function(dds_results, dds_rna, cond_data, exclList=c(), prior=F){
  RNAsp_dds <- list()
  
  
  message("updating dispersions for:")
  message('RNA')
  dispersions(dds_rna[['RNA']])[which(is.na(dispersions(dds_rna[['RNA']])))] <- 10 #max(dispersions(dds_results))
  mcols(dds_results)$dispersion <- dispersions(dds_rna[['RNA']])
  dds_results <- nbinomWaldTest(dds_results, betaPrior = prior)
  RNAsp_dds[['RNA']] <- dds_results
  
  message(paste(dim(dds_results), collapse = "\t"))
  return(RNAsp_dds)
}

expandDups <- function(output){
  output_orig <- output
  if(class(output_orig)=="matrix"){
    output_orig <- as.data.frame(output_orig)
  }
  output_new <- cbind(rownames(output_orig), output_orig)
  colnames(output_new)[1] <- "Row.names"
  # Identify duplicates, if any exist
  message("identifying duplicates")
  dups <- output_new[grep("\\(.*\\)$",output_new$Row.names),]
  dups$Row.names <- gsub("^\\((.*)\\)$","\\1",dups$Row.names)
  # Add everything but the duplicates to the final output
  message("resolving duplicates")
  output_final <- output_new[-(grep("\\(.*\\)$",output_new$Row.names)),]
  output_final <- output_final[!(is.na(output_final$Row.names)),]
  if(nrow(dups) > 0) {
    for(i in 1:nrow(dups)){
      dup_id <- unlist(strsplit(dups$Row.names[i],";"))
      dup_exp <- dups[rep(i, length(dup_id)), ]
      dup_exp$Row.names <- dup_id
      output_final <- rbind(dup_exp,output_final)
    }
    rownames(output_final) <- output_final[,1]
    output_final <- output_final[,-1]
  }else {
    output_final <- output_orig
  }
  return(output_final)
}

cellSpecificTtest<-function(attributesData, counts_norm, dups_output, ctrl_mean, exp_mean, ctrl_cols, exp_cols, altRef = T, correction="BH", cutoff=0.01, prior=F){
  snp_data <- subset(attributesData,allele=="ref" | allele=="alt")
  
  cat("snp_data_In_function_0\n")
  cat(str(snp_data))
  cat("\n")
  
  snp_data$comb <- paste(snp_data$SNP,"_",snp_data$window,"_",snp_data$strand,sep="")
  tmp_ct <- as.data.frame(table(snp_data$comb))
  
  cat("tmp_ct_In_function_0\n")
  cat(str(tmp_ct))
  cat("\n")
  
  snp_data_pairs <- snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq==2,]$Var1,]
  
  cat("snp_data_pairs_In_function_0\n")
  cat(str(snp_data_pairs))
  cat("\n")
  
  snp_data_rejected <- snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq!=2,]$Var1,]
  
  cat("snp_data_rejected_In_function_0\n")
  cat(str(snp_data_rejected))
  cat("\n")
  
  snp_data_ctdata_pairs <- merge(snp_data_pairs, counts_norm, by.x="ID", by.y="row.names", all.x=T, no.dups=F)
  
  cat("snp_data_ctdata_pairs_In_function_0\n")
  cat(str(snp_data_ctdata_pairs))
  cat("\n")
  
  # snp_data_ctdata_pairs$haplotype
  
  snp_data_ctdata_pairs <- snp_data_ctdata_pairs[order(snp_data_ctdata_pairs$SNP, snp_data_ctdata_pairs$window, snp_data_ctdata_pairs$strand, snp_data_ctdata_pairs$allele),]
  
  cat("snp_data_pairs_In_function_1\n")
  cat(str(snp_data_pairs))
  cat("\n")
  
  
  snp_data_expdata_pairs <- merge(snp_data_pairs,dups_output, by.x="ID", by.y="row.names", all.x=T, no.dups=F)
  
  cat("snp_data_expdata_pairs_In_function_0\n")
  cat(str(snp_data_expdata_pairs))
  cat("\n")
  
  snp_data_expdata_pairs <- snp_data_expdata_pairs[order(snp_data_expdata_pairs$SNP, snp_data_expdata_pairs$window, snp_data_expdata_pairs$strand, snp_data_expdata_pairs$allele),]
  
  cat("snp_data_expdata_pairs_In_function_1\n")
  cat(str(snp_data_expdata_pairs))
  cat("\n")
  
  cat("snp_data_pairs_In_function_REMEMBER\n")
  cat(str(snp_data_pairs))
  cat("\n")
  
  if(altRef==T){
    
    cat("Hello_world_1\n")
    
    evens <- seq(2, nrow(snp_data_pairs), by = 2)
    odds <- seq(1, nrow(snp_data_pairs), by = 2)
  }else{
    
    cat("Hello_world_2\n")
    
    evens <- seq(1, nrow(snp_data_pairs), by=2)
    odds <- seq(2, nrow(snp_data_pairs), by=2)
  }
  
  cat("evens_0\n")
  cat(str(evens))
  cat("\n")
  
  cat("odds_0\n")
  cat(str(odds))
  cat("\n")
  
  out <- snp_data_pairs[which(snp_data_pairs$allele=="ref"),c("ID","SNP","chr","pos","ref","alt","allele","window","strand")]
  
  cat("out_In_function_1\n")
  cat(str(out))
  cat("\n")
  
  out$comb <- paste0(out$SNP,"_",out$window,"_",out$strand)
  
  cat("out_In_function_2\n")
  cat(str(out))
  cat("\n")
  
  cat("snp_data_expdata_pairs_In_function_REMEMBER\n")
  cat(str(snp_data_expdata_pairs))
  cat("\n")
  
  
  merge_companion_REF<-snp_data_expdata_pairs[which(snp_data_expdata_pairs$allele=="ref"),c("ID","comb","ctrl_mean","exp_mean","log2FoldChange","lfcSE","pvalue","padj")]
  
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'ctrl_mean')]<-"A_Ctrl_Mean"
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'exp_mean')]<-"A_Exp_Mean"
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'log2FoldChange')]<-"A_log2FC"
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'lfcSE')]<-"A_log2FC_SE"
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'pvalue')]<-"A_logP"
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'padj')]<-"A_logPadj_BH"
  
  cat("merge_companion_REF_In_function_0\n")
  cat(str(merge_companion_REF))
  cat("\n")
  
  
  out <- merge(out, merge_companion_REF, by = c("comb","ID"))
  
  cat("out_In_function_3_POST_merge_with_merge_companion\n")
  cat(str(out))
  cat("\n")
 
  
  out$A_logPadj_BH <- -log10(out$A_logPadj_BH)
  out$A_logPadj_BH[out$A_logPadj_BH < 0] <- 0
  out$A_logPadj_BH[out$A_logPadj_BH == Inf] <- max(out$A_logPadj_BH[is.finite(out$A_logPadj_BH)])
  out$A_logPadj_BF <- -log10(out$A_logP*(nrow(snp_data_expdata_pairs)/2))
  out$A_logP <- -log10(out$A_logP)
  out$A_logP[is.na(out$A_logP)] <- 0
  out$A_logP[out$A_logP == Inf] <- max(out$A_logP[is.finite(out$A_logP)])
  out$A_logPadj_BF[out$A_logPadj_BF < 0] <- 0
  out$A_logPadj_BF[out$A_logPadj_BF == Inf] <- max(out$A_logPadj_BF[is.finite(out$A_logPadj_BF)])
  
  merge_companion_ALT<-snp_data_expdata_pairs[which(snp_data_expdata_pairs$allele=="alt"),c("comb","ctrl_mean","exp_mean","log2FoldChange","lfcSE","pvalue","padj")]
  
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'ctrl_mean')]<-"B_Ctrl_Mean"
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'exp_mean')]<-"B_Exp_Mean"
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'log2FoldChange')]<-"B_log2FC"
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'lfcSE')]<-"B_log2FC_SE"
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'pvalue')]<-"B_logP"
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'padj')]<-"B_logPadj_BH"
  
  
  cat("merge_companion_ALT_In_function_0\n")
  cat(str(merge_companion_ALT))
  cat("\n")
  
  
  out <- merge(out, merge_companion_ALT, by = "comb")
  
 

  cat("out_In_function_6\n")
  cat(str(out))
  cat("\n")
  
  out$B_logPadj_BH <- -log10(out$B_logPadj_BH)
  out$B_logPadj_BH[out$B_logPadj_BH < 0] <- 0
  out$B_logPadj_BH[out$B_logPadj_BH == Inf] <- max(out$B_logPadj_BH[is.finite(out$B_logPadj_BH)])
  out$B_logPadj_BF <- -log10(out$B_logP*(nrow(snp_data_expdata_pairs)/2))
  out$B_logP <- -log10(out$B_logP)
  out$B_logP[is.na(out$B_logP)] <- 0
  out$B_logP[out$B_logP == Inf] <- max(out$B_logP[is.finite(out$B_logP)])
  out$B_logPadj_BF[out$B_logPadj_BF < 0] <- 0
  out$B_logPadj_BF[out$B_logPadj_BF == Inf] <- max(out$B_logPadj_BF[is.finite(out$B_logPadj_BF)])
  
  cat("out_In_function_7\n")
  cat(str(out))
  cat("\n")
  
  out2 <- out#[,c(1:12, 16:19, 26:28, 23:25)]
  
  cat("out2_in_function_0\n")
  cat(str(out2))
  cat("\n")
  
  
  
  # Don't try to do the t test for ones with all zeros.
  ignore_idx <- which(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[odds, exp_cols]) < 10 |
                        is.na(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]))  | is.na(rowMeans(snp_data_ctdata_pairs[odds, exp_cols])) |
                        rowMeans(snp_data_ctdata_pairs[evens,ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[evens, exp_cols]) < 10 |
                        is.na(rowMeans(snp_data_ctdata_pairs[evens,ctrl_cols]))  | is.na(rowMeans(snp_data_ctdata_pairs[evens, exp_cols])) )
  
  # For the numerator, set zero values to 1 so that the log-ratio is defined.
  counts1 <- snp_data_ctdata_pairs
  counts1[counts1 == 0] <- 1
  
  # t test
  ratios_A <- log2((counts1[evens, exp_cols]) / rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols]))
  ratios_B <- log2((counts1[odds, exp_cols]) / rowMeans(snp_data_ctdata_pairs[odds, ctrl_cols]))
  
  ratios_list <- list(ratios_A, ratios_B)
  
  t_pvalue <- sapply(1:nrow(ratios_A), function(i) if(i %in% ignore_idx){NA} else{
    t.test(as.numeric(ratios_A[i,]), as.numeric(ratios_B[i,]), var.equal=F, paired=T)$p.value})
  t_stat <- sapply(1:nrow(ratios_A), function(i) if(i %in% ignore_idx){NA} else{
    t.test(as.numeric(ratios_A[i,]), as.numeric(ratios_B[i,]), var.equal=F, paired=T)$statistic})
  
  if(prior==T){
    mean_ratios_A <- log2(rowMeans(counts1[evens,exp_cols], na.rm = T) / rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols], na.rm = T))
    mean_ratios_B <- log2(rowMeans(counts1[odds,exp_cols], na.rm = T) / rowMeans(snp_data_ctdata_pairs[odds, ctrl_cols], na.rm = T))
    
    out2$Log2Skew <- mean_ratios_B - mean_ratios_A
  }
  
  if(prior==F){
    out2$Log2Skew <- out2$B_log2FC - out2$A_log2FC
  }
  out2$LogSkew_SE <- sqrt(out2$A_log2FC_SE^2+out2$B_log2FC_SE^2)
  out2$Skew_logP <- ifelse(is.na(t_pvalue), 0, -log10(t_pvalue))
  
  OE_threshold <- -log10(cutoff)
  if(correction=="BF"){
    is_OE <- out2$A_logPadj_BF >= OE_threshold | out2$B_logPadj_BF >= OE_threshold
  }
  if(correction=="BH"){
    is_OE <- out2$A_logPadj_BH >= OE_threshold | out2$B_logPadj_BH >= OE_threshold
  }
  out2$Skew_logFDR <- -log10(p.adjust(t_pvalue, method = "BH"))
  out2$Skew_logFDR_act <- rep(NA, nrow(out))
  q_idx <- intersect(which(is_OE), which(!is.na(t_pvalue)))
  out2$Skew_logFDR_act[q_idx] <- -log10(p.adjust(t_pvalue[q_idx],method="BH"))
  
  cat("out2_in_function_1\n")
  cat(str(out2))
  cat("\n")
  
  return(out2)
}

DESkew <- function(conditionData, counts_norm, attributesData, celltype, dups_output,prior, cutoff, paired=F){
  
  ds_cond_data <- as.data.frame(conditionData[which(conditionData$condition=="DNA" | conditionData$condition==celltype),,drop=F])
  
  cat("ds_cond_data_in_function_Post_merge_merge_companion_REF\n")
  cat(str(ds_cond_data))
  cat("\n")
  
  # Prepare the sample table
  message("Preparing Sample Table")
  dna_reps <- nrow(as.data.frame(ds_cond_data[which(ds_cond_data$condition=="DNA"),]))
  rna_reps <- nrow(as.data.frame(ds_cond_data[which(ds_cond_data$condition==celltype),]))
  avg_reps <- (dna_reps+rna_reps)/2
  total_cond <- length(unique(ds_cond_data$condition))
  samps <- data.frame(material=factor(c(rep("DNA",dna_reps*total_cond),rep("RNA",rna_reps*total_cond))),
                      allele=factor(rep(c("ref","alt"),((dna_reps+rna_reps))), levels = c("ref","alt")),
                      sample=factor(rep(c(rownames(ds_cond_data)[which(ds_cond_data$condition=="DNA")], rownames(ds_cond_data)[which(ds_cond_data$condition==celltype)]),each=total_cond)))
  
  cat("samps_in_function_Post_merge_merge_companion_REF\n")
  cat(str(samps))
  cat("\n")
  
  # Reorganize the count data
  message("reorganizing count data")
  snp_data <- subset(attributesData,allele=="ref" | allele=="alt")
  
  cat("snp_data_in_function_Post_merge_merge_companion_REF\n")
  cat(str(snp_data))
  cat("\n")
  
  snp_data$comb <- paste(snp_data$SNP,"_",snp_data$window,"_",snp_data$strand,sep="")
  
  tmp_ct <- as.data.frame(table(snp_data$comb))
  
  snp_data_pairs <- snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq==2,]$Var1,]
  
  cat("snp_data_pairs_in_function_Post_merge_merge_companion_REF\n")
  cat(str(snp_data_pairs))
  cat("\n")
  
  pairs_ids <- snp_data_pairs$ID
  
  cat("dups_output_in_function_REMEMBER\n")
  cat(str(dups_output))
  cat("\n")
  
  dups_ids <- rownames(dups_output)
  
  snp_data_pairs <- merge(snp_data_pairs,
                          dups_output, by.x="ID", by.y="row.names", all=T, no.dups=F)
  
  cat("snp_data_pairs_in_function_Post_merge_dups_output\n")
  cat(str(snp_data_pairs))
  cat("\n")
  
  message("COUNT OUT")
  
  out <- snp_data_pairs[which(snp_data_pairs$allele=="ref"),c("ID","SNP","chr","pos","ref","alt","allele","window","strand")]
  
  cat("out_in_function_0\n")
  cat(str(out))
  cat("\n")
  
  out$comb <- paste0(out$SNP,"_",out$window,"_",out$strand,sep="")
  
 
  
  merge_companion_REF<-snp_data_pairs[which(snp_data_pairs$allele=="ref"),c("ID","comb","ctrl_mean","exp_mean","log2FoldChange","lfcSE","pvalue","padj")]
  
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'ctrl_mean')]<-"A_Ctrl_Mean"
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'exp_mean')]<-"A_Exp_Mean"
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'log2FoldChange')]<-"A_log2FC"
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'lfcSE')]<-"A_log2FC_SE"
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'pvalue')]<-"A_logP"
  colnames(merge_companion_REF)[which(colnames(merge_companion_REF) == 'padj')]<-"A_logPadj_BH"
  
  cat("merge_companion_REF_In_function_0\n")
  cat(str(merge_companion_REF))
  cat("\n")
  
  
  out <- merge(out, merge_companion_REF, by = c("comb","ID"))
  
  cat("out_in_function_Post_merge_merge_companion_REF\n")
  cat(str(out))
  cat("\n")

  out$A_logPadj_BH <- -log10(out$A_logPadj_BH)
  out$A_logPadj_BH[out$A_logPadj_BH < 0] <- 0
  out$A_logPadj_BH[out$A_logPadj_BH == Inf] <- max(out$A_logPadj_BH[is.finite(out$A_logPadj_BH)])
  out$A_logPadj_BF <- -log10(out$A_logP*(nrow(snp_data_pairs)/2))
  out$A_logP <- -log10(out$A_logP)
  out$A_logP[is.na(out$A_logP)] <- 0
  out$A_logP[out$A_logP == Inf] <- max(out$A_logP[is.finite(out$A_logP)])
  out$A_logPadj_BF[out$A_logPadj_BF < 0] <- 0
  out$A_logPadj_BF[out$A_logPadj_BF == Inf] <- max(out$A_logPadj_BF[is.finite(out$A_logPadj_BF)])
  
  merge_companion_ALT<-snp_data_pairs[which(snp_data_pairs$allele=="alt"),c("comb","ctrl_mean","exp_mean","log2FoldChange","lfcSE","pvalue","padj")]
  
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'ctrl_mean')]<-"B_Ctrl_Mean"
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'exp_mean')]<-"B_Exp_Mean"
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'log2FoldChange')]<-"B_log2FC"
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'lfcSE')]<-"B_log2FC_SE"
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'pvalue')]<-"B_logP"
  colnames(merge_companion_ALT)[which(colnames(merge_companion_ALT) == 'padj')]<-"B_logPadj_BH"
  
  cat("merge_companion_ALT_In_function_0\n")
  cat(str(merge_companion_ALT))
  cat("\n")
  
  
  out <- merge(out, merge_companion_ALT, by = "comb")
  
  cat("out_in_function_Post_merge_merge_companion_ALT\n")
  cat(str(out))
  cat("\n")

  out$B_logPadj_BH <- -log10(out$B_logPadj_BH)
  out$B_logPadj_BH[out$B_logPadj_BH < 0] <- 0
  out$B_logPadj_BH[out$B_logPadj_BH == Inf] <- max(out$B_logPadj_BH[is.finite(out$B_logPadj_BH)])
  out$B_logPadj_BF <- -log10(out$B_logP*(nrow(snp_data_pairs)/2))
  out$B_logP <- -log10(out$B_logP)
  out$B_logP[is.na(out$B_logP)] <- 0
  out$B_logP[out$B_logP == Inf] <- max(out$B_logP[is.finite(out$B_logP)])
  out$B_logPadj_BF[out$B_logPadj_BF < 0] <- 0
  out$B_logPadj_BF[out$B_logPadj_BF == Inf] <- max(out$B_logPadj_BF[is.finite(out$B_logPadj_BF)])
  
  counts_allele_id <- merge(snp_data_pairs[,c("ID","comb","allele","SNP","window","strand")], counts_norm, by.x="ID", by.y="row.names", all.x=T)
  
  rownames(counts_allele_id) <- counts_allele_id$ID
  
  cat("counts_allele_id_in_function_0\n")
  cat(str(counts_allele_id))
  cat("\n")
  
  message("done")
  
  message("COUNTS REF")
  
  counts_ref_all <- counts_allele_id[which(counts_allele_id$allele=="ref"),,drop=F]
  
  cat("counts_ref_all_in_function_0\n")
  cat(str(counts_ref_all))
  cat("\n")
  
  counts_ref <- counts_ref_all[,colnames(counts_ref_all) %in% rownames(ds_cond_data),drop=F]
  
  cat("counts_ref_in_function_0\n")
  cat(str(counts_ref))
  cat("\n")
  
  colnames(counts_ref) <- paste0(colnames(counts_ref),"_ref")
  counts_ref$comb <- paste0(counts_ref_all$SNP,"_",counts_ref_all$window,"_",counts_ref_all$strand)
  counts_ref$ID <- rownames(counts_ref)
  
  counts_ref <- merge(counts_ref, snp_data_pairs[which(snp_data_pairs$ID %in% rownames(counts_ref)),c("ID","SNP","chr","pos","ref","alt","allele","window","strand","comb")],by=c("ID","comb"),all.x=T)
  
  cat("counts_ref_in_function_1\n")
  cat(str(counts_ref))
  cat("\n")
  
  message("done")  
  
  message("COUNTS ALT")
  
  counts_alt_all <- counts_allele_id[which(counts_allele_id$allele=="alt"),,drop=F]
  
  cat("counts_alt_all_in_function_0\n")
  cat(str(counts_alt_all))
  cat("\n")
  
  counts_alt <- counts_alt_all[,colnames(counts_alt_all) %in% rownames(ds_cond_data),drop=F]
  
  cat("counts_alt_in_function_0\n")
  cat(str(counts_alt))
  cat("\n")
  
  colnames(counts_alt) <- paste0(colnames(counts_alt),"_alt")
  counts_alt$comb <- paste0(counts_alt_all$SNP,"_",counts_alt_all$window,"_",counts_alt_all$strand)
  counts_alt$ID <- rownames(counts_alt)
  counts_alt <- merge(counts_alt, snp_data_pairs[which(snp_data_pairs$ID %in% rownames(counts_alt)),c("ID","SNP","chr","pos","ref","alt","allele","window","strand","comb")],by=c("ID","comb"),all.x=T)
  
  cat("counts_alt_in_function_1\n")
  cat(str(counts_alt))
  cat("\n")
  
  message("done")
  # message(paste0(colnames(counts_alt), collapse = "\t"))
  # message(paste0(colnames(counts_ref), collapse = "\t"))
  
  counts_ref_alt <- merge(counts_ref, counts_alt, by=c("SNP","chr","pos","ref","alt","window","strand","comb"), all=T)
  
  cat("counts_ref_alt_in_function_POST_merge_counts_ref_and_counts_alt\n")
  cat(str(counts_ref_alt))
  cat("\n")
  
  message("counts merged") 
  
  column_order <- data.frame(allele=factor(rep(c("ref","alt"),((dna_reps+rna_reps))), levels = c("ref","alt")),
                             sample=factor(rep(c(rownames(ds_cond_data)[which(ds_cond_data$condition=="DNA")], rownames(ds_cond_data)[which(ds_cond_data$condition==celltype)]),each=total_cond)))
  column_order$order <- paste0(column_order$sample,"_",column_order$allele)
  
  cat("column_order_in_function_0\n")
  cat(str(column_order))
  cat("\n")
  
  message("order set")
  message(paste0(colnames(counts_ref_alt), collapse = "\t"))
  counts_ref_alt <- counts_ref_alt[,c("ID.x","SNP","chr","pos","ref","alt","allele.x","window","strand",column_order$order)]
  
  cat("counts_ref_alt_in_function_1\n")
  cat(str(counts_ref_alt))
  cat("\n")
  
  colnames(counts_ref_alt) <- c("ID","SNP","chr","pos","ref","alt","allele","window","strand",column_order$order)
  message(paste0(dim(counts_ref_alt), collapse = "\t"))
  
  cat("counts_ref_alt_in_function_2\n")
  cat(str(counts_ref_alt))
  cat("\n")
  
  counts_comp <- counts_ref_alt[complete.cases(counts_ref_alt[,column_order$order]),]
  
  cat("counts_comp_in_function_0\n")
  cat(str(counts_comp))
  cat("\n")
  
  counts_mat <- as.matrix(counts_comp[,column_order$order])
  
  cat("counts_mat_in_function_0\n")
  cat(str(counts_mat))
  cat("\n")
  
  ids_comp <- counts_comp$ID
  
  cat("ids_comp_in_function_0\n")
  cat(str(ids_comp))
  cat("\n")
  
  message("ids pulled")
  
  # Set Design definition
  design <- ~sample + allele
  
  # Run DESeq analysis
  message("Running DESeq")
  
  cat("counts_mat_in_function_0_to_run_DESeqDataSetFromMatrix\n")
  cat(str(counts_mat))
  cat("\n")
  
  cat("samps\n")
  cat(str(samps))
  cat("\n")
  
  cat("design\n")
  cat(str(design))
  cat("\n")
  
  dds <- DESeqDataSetFromMatrix(counts_mat, samps, design)
  
  cat("dds_in_function_0\n")
  cat(str(dds))
  cat("\n")
  
  if(paired==F){
    design(dds) <- ~material + allele + material:allele
    
    dds <- DESeq(dds, fitType = "local", minReplicatesForReplace = Inf)
    
    res.diff <- results(dds, name="materialRNA.allelealt",cooksCutoff=F,independentFiltering=F)
    
  }
  
  if(paired==T){
    message("running paired samples")
    
    sample_lets <- c(rep(LETTERS[1:dna_reps], each=total_cond), rep(LETTERS[1:rna_reps], each=total_cond))
    dds$sample.n <- as.factor(sample_lets)
    
    design(dds) <- ~material + material:sample.n + material:allele
    sizeFactors(dds) <- rep(1, (dna_reps+rna_reps)*total_cond)
    
    # if(dna_reps != rna_reps){
    #   warning("Number of DNA replicates and RNA replicates unequal. Using the lower number of replicates to run paired samples. If you want to avoid this run MPRAmodel with `paired=F`")
    #   mm <- model.matrix(~material + material:sample.n + material:allele, colData(dds))
    #   col_mm <- ncol(mm)
    #   mm <- mm[,c(1:((min(dna_reps,rna_reps))*2),(col_mm-1),col_mm)]
    #   dds <- DESeq(dds, full = mm, fitType = "local", minReplicatesForReplace=Inf)
    # }
    # else{
    dds <- DESeq(dds, fitType = "local", minReplicatesForReplace = Inf)
    # }
    
    #Get the skew results
    # cell_res <- paste0("condition",celltype,".countalt")
    message(paste0(resultsNames(dds), collapse = "\t"))
    
    res.diff <- results(dds, contrast=list("materialRNA.allelealt","materialDNA.allelealt"), cooksCutoff=F, independentFiltering=F)
  }
  
  # res.diff <- results(dds, name="materialRNA.allelealt",cooksCutoff=F,independentFiltering=F)
  
  res.diff <- as.data.frame(res.diff)[,-1]
  
  cat("res.diff_in_function_0\n")
  cat(str(res.diff))
  cat("\n")
  
  colnames(res.diff) <- c("Log2Skew","Skew_SE","skewStat","Skew_logP","Skew_logFDR")
  
  cat("res.diff_in_function_1\n")
  cat(str(res.diff))
  cat("\n")
  
  res.diff$Skew_logP <- -log10(as.data.frame(res.diff)$Skew_logP)
  res.diff$Skew_logFDR <- -log10(as.data.frame(res.diff)$Skew_logFDR)
  message(paste0("res_diff samples: ", nrow(res.diff)))
  res.diff$ID <- ids_comp
  
  cat("res.diff_in_function_2\n")
  cat(str(res.diff))
  cat("\n")
  
  message("combining data")
  res_comp <- merge(out, res.diff, by="ID", all.y=T)
  
  cat("res_comp_in_function_0\n")
  cat(str(res_comp))
  cat("\n")
  
  OE_threshold <- -log10(cutoff)
  
  cat("OE_threshold_in_function_0\n")
  cat(str(OE_threshold))
  cat("\n")
  
  glm_stats <- res_comp
  
  cat("glm_stats_in_function_0\n")
  cat(str(glm_stats))
  cat("\n")
  
  glm_stats$p_orig <- 10^(-glm_stats$Skew_logP)
  glm_stats$OE <- ifelse(glm_stats$A_logPadj_BH > OE_threshold | glm_stats$B_logPadj_BH > OE_threshold, T, F)
  
  cat("glm_stats_in_function_1\n")
  cat(str(glm_stats))
  cat("\n")
  
  glm_update <- glm_stats[which(glm_stats$OE),c("ID","p_orig")]
  glm_update$Skew_logFDR_act <- -log10(p.adjust(glm_update$p_orig, method = "BH"))
  
  cat("glm_stats_in_function_2\n")
  cat(str(glm_stats))
  cat("\n")
  
  res_comp <- merge(res_comp, glm_update[,c("ID","Skew_logFDR_act")], by="ID", all=T)
  
  cat("res_comp_in_function_FINAL\n")
  cat(str(res_comp))
  cat("\n")
  
  return(res_comp)
}

### Retrieve output data for future functions - if only looking for results and not the plots this is the only function that needs to be called
# Any subsequent functions should only need an output from here
## INPUT
# countsData      : table of tag counts, columns should include: Barcode, Oligo, Sample names
# attributesData  : table of full attributes, columns should include: ID, SNP, Project, Window, Strand, Allele,
# Haplotype, Bash. **NB** If you are just running this function make sure to pass your attributes table through the addHaplo function
# conditionData   : table of conditions, 2 columns no header align to the variable column headers of countsData
# and 'RNA'
# exclList        : list of celltypes that should be excluded from the analysis. Empty by default
# altRef          : LOGICAL default T, indicating the order to sort alleles for allelic skew. alt/ref default
# file_prefix     : String to indicate what written file names include
# method          : Method to use for normalization
# negCtrlName     : String indicating what negative controls are called in the attributes table
# tTest           : LOGICAL default T, identify emVARs using the tTest method
# DEase           : LOGICAL default T, identify emVARS using the DESeq method of determining allelic skew
# correction      : String indicating whether to use Benjamini Hochberg ("BH", default) or Bonferroni ("BF") for p-value correction
# cutoff          : significance cutoff for including alleles for skew calculation (tTest only)
# upDisp          : LOGICAL default T, update dispersions with 'RNA' specific calculations
# prior           : LOGICAL default T, use betaPrior=T when calculating the 'RNA' specific dispersions.
## OUTPUT: writes duplicate output and ttest files for each 'RNA'

# Bug order attributesData, conditionData, -> conditionData, attributesData

dataOut <- function(countsData, conditionData,attributesData, exclList = c(), altRef = T, file_prefix, method = 'ss',negCtrlName="negCtrl",
                    tTest=T, DEase=T, cSkew=T, correction="BH", cutoff=0.01, upDisp=T, prior=F, paired=F){
  
  message("------------------------------------------------------------------------------->starting_datasets")
  

  cat("countsData_In_function_0\n")
  cat(str(countsData))
  cat("\n")
  
  cat("attributesData_In_function_0\n")
  cat(str(attributesData))
  cat("\n")
  
  cat("conditionData_In_function_0\n")
  cat(str(conditionData))
  cat("\n")
  
  cat("file_prefix_In_function_0\n")
  cat(sprintf(as.character(file_prefix)))
  cat("\n")

  
  dds_results_all <- tagNorm(countsData, conditionData, attributesData, exclList, method, negCtrlName, upDisp, prior)
  
  # cat("dds_results_all_In_function_0\n")
  # cat(str(dds_results_all))
  # cat("\n")
  
  dds_results <- dds_results_all[[1]]
  dds_results_orig <- dds_results_all[[2]]
  attribute_ids <- (attributesData[attributesData$ctrl_exp==negCtrlName,])$ID
  
  cat("attribute_ids_In_function_0\n")
  cat(str(attribute_ids))
  cat("\n")
  
  full_attribute_ids <- attributesData$ID
  
  cat("full_attribute_ids_In_function_0\n")
  cat(str(full_attribute_ids))
  cat("\n")
  
  message("Tags Normalized")
  cond_data <- conditionStandard(conditionData)
  
  cat("cond_data_In_function_0\n")
  cat(str(cond_data))
  cat("\n")
  cat(sprintf(as.character(row.names(cond_data))))
  cat("\n")
  
  message("condition data standardized")
  
  count_data <- oligoIsolate(countsData)
  
  cat("count_data_In_function_Second time\n")
  cat(str(count_data))
  cat("\n")
  
  message("Oligos isolated Second time")
  
  
  colnames(count_data) <- row.names(cond_data)
  
  full_output<-list()
  full_output_var<-list()
  
  condition_table <- as.data.frame(cond_data)
  
  cat("condition_table_In_function_0\n")
  cat(str(condition_table))
  cat("\n")
  
  dds_rna <- list()
  # Celltype based DESeq analysis
  
  rna_cols <- cond_data[which(cond_data$condition=='RNA'),]
  rna_count <- count_data[,rownames(rna_cols)]
  dds_rna_temp <- DESeqDataSetFromMatrix(countData = rna_count, colData = rna_cols, design = ~1)
  sizeFactors(dds_rna_temp) <- sizeFactors(dds_results)[rownames(rna_cols)]
  dds_rna_temp <- estimateDispersions(dds_rna_temp, fitType='local')
  dds_rna[['RNA']] <- dds_rna_temp
  
  # cat("dds_rna_temp_In_function_0\n")
  # cat(str(dds_rna_temp))
  # cat("\n")
  
  

  
  # Replace dispersions in normalized dds with the 'RNA' specific dispersions
  if(upDisp==T){
    RNAsp_dds <- tagSig(dds_results, dds_rna, cond_data, exclList, prior)
  }
  
  # cat("RNAsp_dds_In_function_0\n")
  # cat(str(RNAsp_dds))
  # cat("\n")
  
  if(upDisp==T){
    dds_results <- RNAsp_dds[['RNA']]
  }
  
  #### Print normalized counts
  
  counts_norm_all <- counts(dds_results, normalized = T)
  
  cat("counts_norm_all_In_function_0\n")
  cat(str(counts_norm_all))
  cat("\n")
  
  # write.table(counts_norm_all,paste0("results/", file_prefix, "_", fileDate(),"_normalized_counts.out"), quote = F, sep = "\t")
  write.table(counts_norm_all,paste0("results/", file_prefix,"_normalized_counts.out"), quote = F, sep = "\t")
  
  
  counts_norm <- counts(dds_results, normalized = T)
  counts_norm <- counts_norm[,colnames(counts_norm) %in% rownames(cond_data)[which(condition_table$condition %in% c("DNA",'RNA'))]]
  
  cat("counts_norm_RNA????In_function_0\n")
  cat(str(counts_norm))
  cat("\n")
  
  # write.table(counts_norm,paste0("results/", file_prefix, "_", fileDate(),"_",'RNA', "_normalized_counts.out"), quote = F, sep = "\t")
  
  write.table(counts_norm,paste0("results/", file_prefix,"_",'RNA', "_normalized_counts.out"), quote = F, sep = "\t")
  
  if(DEase==T){
    message("Removing count duplicates")
    counts_DE <- counts(dds_results, normalized=F)
    counts_norm_DE <- expandDups(counts_DE)
  }
  
  outputA <- results(dds_results, contrast=c("condition",'RNA',"DNA"), cooksCutoff=F, independentFiltering=F)
  
  cat("outputA_In_function_0\n")
  cat(str(outputA))
  cat("\n")
  
  temp_outputB <- results(dds_results_orig, contrast=c("condition",'RNA',"DNA"), cooksCutoff=F, independentFiltering=F)
  
  cat("temp_outputB_In_function_0\n")
  cat(str(temp_outputB))
  cat("\n")
  
  # message("Plotting Normalization Curves")
  # # pdf(paste0("plots/",file_prefix,"_",fileDate(),"_Normalized_FC_Density_",'RNA',".pdf"),width=10,height=10)
  # pdf(paste0("plots/",file_prefix,"_Normalized_FC_Density_",'RNA',".pdf"),width=10,height=10)
  # 
  # plot(density(temp_outputB[attribute_ids,"log2FoldChange"],na.rm=TRUE),xlim=c(-3,3),ylim=c(0,1.5),col="grey",main=paste0("Normalization - ",'RNA'))
  # lines(density(temp_outputB$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="black")
  # lines(density(outputA$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="red")
  # lines(density(outputA[attribute_ids,"log2FoldChange"],na.rm=TRUE),xlim=c(-3,3),col="salmon")
  # text(1.5,0.4,adj=c(0,0),labels="All - baseline",col="black")
  # text(1.5,0.35,adj=c(0,0),labels="All - corrected",col="red")
  # text(1.5,0.3,adj=c(0,0),labels=paste0(negCtrlName," - baseline"),col="grey")
  # text(1.5,0.25,adj=c(0,0),labels=paste0(negCtrlName," - corrected"),col="salmon")
  # abline(v=0)
  # dev.off()
  
  message("Results of dds_results recieved")
  
  ctrl_cols <- row.names(condition_table[condition_table$condition=="DNA",])
  exp_cols <- row.names(condition_table[condition_table$condition=='RNA',])
  
  message("Control and Experiment Columns Set")
  
  DNA_mean <- rowMeans(count_data[, colnames(count_data) %in% ctrl_cols], na.rm=T)
  ctrl_mean <- rowMeans(counts_norm[, colnames(counts_norm) %in% ctrl_cols], na.rm = T)
  exp_mean <- rowMeans(counts_norm[, colnames(counts_norm) %in% exp_cols], na.rm = T)
  output_2 <- cbind(DNA_mean,ctrl_mean,exp_mean,outputA[,-1])
  
  counts_norm_sp <- expandDups(counts_norm)
  
  message("Removing Duplicates")
  dups_output<-expandDups(output_2)
  # message(paste0(colnames(dups_output),collapse = "\t"))
  
  message("Writing Standard Results File")
  full_outputA<-merge(attributesData, as.data.frame(dups_output), by.x="ID", by.y="row.names", all.x=TRUE, no.dups=F)
  full_output[['RNA']]<-full_outputA
  # write.table(full_outputA, paste0("results/", file_prefix, "_", 'RNA', "_", fileDate(), ".out"), row.names=F, col.names=T, sep="\t", quote=F)
  write.table(full_outputA, paste0("results/", file_prefix, "_", 'RNA',".out"), row.names=F, col.names=T, sep="\t", quote=F)
  
  # ############################################################################################################################################################
  # message("----------------------------------------------------------------------------------------------------------->Next is cellSpecificTtest")
  # ############################################################################################################################################################
  # 
  # 
  # 
  # if(tTest==T){
  #   message("Writing T-Test Results File")
  #   outA<-cellSpecificTtest(attributesData, counts_norm_sp, dups_output, ctrl_mean, exp_mean, ctrl_cols, exp_cols, altRef, correction, cutoff, prior)
  #   
  #   cat("outA_In_function_0\n")
  #   cat(str(outA))
  #   cat("\n")
  #   
  #   
  #   full_output_var[['RNA']]<-outA
  #   # write.table(outA,paste0("results/", file_prefix, "_", 'RNA', "_emVAR_ttest_", fileDate(),".out"), row.names=F, col.names=T, sep="\t", quote=F)
  #   write.table(outA,paste0("results/", file_prefix, "_", 'RNA', "_emVAR_ttest", ".out"), row.names=F, col.names=T, sep="\t", quote=F)
  #   # write.table(attributesData, "results/ttest_attributes.txt", row.names = F, quote = F, sep = "\t")
  #   # write.table(counts_norm_sp, paste0("results/ttest_",'RNA',"_norm_counts.txt"), quote = F, sep = "\t")
  #   # write.table(dups_output, "results/ttest_dups_out.txt", quote = F, sep = "\t")
  # }
  
   
  
  
  if(DEase==T){
    message("Writing DESeq Allelic Skew Results File")
    if(paired==F){
      cond_pass <- condition_table
    }
    if(paired==T){
      plas_row <- length(conditionData$condition[which(conditionData=="DNA")])
      message(plas_row)
      cell_row <- length(conditionData$condition[which(conditionData$condition=='RNA')])
      message(cell_row)
      if(plas_row==cell_row){
        cond_pass <- condition_table
      }
      else if(plas_row!=cell_row){
        if(plas_row > cell_row){
          drop_num <- plas_row - cell_row
          message(paste0("dropping ", drop_num," DNA replicates for paired GLM analysis"))
          plas_ids <- data.frame(matrix(ncol = 2, nrow=plas_row))
          colnames(plas_ids) <- c("rep","bc_count")
          plas_ids$rep <- rownames(conditionData)[which(conditionData$condition=="DNA")]
          
          for(id in plas_ids$rep){
            # message(id)
            # message(sum(counts_data[,id] > 0))
            
            plas_ids[plas_ids$rep==id,"bc_count"] <- sum(counts_data[,id] > 0)
          }
          plas_ids <- plas_ids[order(plas_ids$bc_count),]
          
          plas_drop <- plas_ids$rep[1:drop_num]
          
          message(paste0("Dropping: ",plas_drop, collapse = "\t"))
          cond_pass <- condition_table[which(!rownames(condition_table)%in%plas_drop),]
        }
        else if(cell_row > plas_row){
          drop_num <- cell_row - plas_row
          message(paste0("dropping ", drop_num," ",'RNA', " replicates for paired GLM analysis"))
          cell_ids <- data.frame(matrix(ncol = 2, nrow=cell_row))
          colnames(cell_ids) <- c("rep","bc_count")
          cell_ids$rep <- rownames(conditionData)[which(conditionData$condition=='RNA')]
          message(paste0(cell_ids$rep, collapse = "\t"))
          
          message(class(cell_ids))
          
          for(id in cell_ids$rep){
            
            cell_ids[cell_ids$rep==id,"bc_count"] <- sum(counts_data[,id] > 0)
          }
          cell_ids <- cell_ids[order(cell_ids$bc_count),]
          
          cell_drop <- cell_ids$rep[1:drop_num]
          
          message(paste0("Dropping: ",cell_drop, collapse = "\t"))
          cond_pass <- condition_table[which(!rownames(condition_table)%in%cell_drop),]
        }
      }
      counts_norm_DE <- counts_norm_DE[,which(colnames(counts_norm_DE) %in% rownames(cond_pass))]
    }
    
    ############################################################################################################################################################
    message("----------------------------------------------------------------------------------------------------------->Next is Allelic Skew")
    ############################################################################################################################################################
    
    
    outB <- DESkew(cond_pass, counts_norm_DE, attributesData, 'RNA', dups_output,prior, cutoff, paired)
    
    cat("-------------------------------------------------------------------------------------------------------------->outB_In_function_0\n")
    cat(str(outB))
    cat("\n")
    
    if(paired==F){
      # write.table(outB,paste0("results/", file_prefix, "_", 'RNA', "_emVAR_glm_", fileDate(),".out"), row.names=F, col.names=T, sep="\t", quote=F)
      write.table(outB,paste0("results/", file_prefix, "_", 'RNA',"_emVAR_glm", ".out"), row.names=F, col.names=T, sep="\t", quote=F)
      write.table(outB,paste(file_prefix, "_", 'RNA',"_emVAR_glm", ".out"), row.names=F, col.names=T, sep="\t", quote=F)
      
    }
    if(paired==T){
      # write.table(outB,paste0("results/", file_prefix, "_", 'RNA', "_emVAR_glm_paired_", fileDate(),".out"), row.names=F, col.names=T, sep="\t", quote=F)
      
      write.table(outB,paste0("results/", file_prefix, "_", 'RNA',"_emVAR_glm_paired",".out"), row.names=F, col.names=T, sep="\t", quote=F)
    }
  }
  
  message("Writing bed File")
  full_bed_outputA<-merge(attributesData, as.matrix(dups_output),by.x="ID",by.y="row.names",all.x=TRUE,no.dups=FALSE)
  
  full_bed_outputA$start <- ""
  full_bed_outputA$stop <- ""
  
  cat("full_bed_outputA_In_function_0\n")
  cat(str(full_bed_outputA))
  cat("\n")
  
  
  # message(paste0(colnames(full_bed_outputA), collapse = "\t"))
  #printbed<-full_bed_outputA[,c("chr","start","stop","ID","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","cigar","md-tag","project")]
  # if(!(c("start","stop") %in% colnames(full_bed_outputA))){
  #   full_bed_outputA$start <- ""
  #   full_bed_outputA$stop <- ""
  # }
  printbed<-full_bed_outputA[,c("chr","start","stop","ID","strand","log2FoldChange","ctrl_mean","exp_mean","pvalue","padj","lfcSE","project")]
  printbed$score<-"."
  
  cat("printbed_In_function_0\n")
  cat(str(printbed))
  cat("\n")
  
  #printbed<-printbed[,c("chr","start","stop","ID","score","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","cigar","md-tag","project")]
  #colnames(printbed)<-c("chr","start","stop","id","score","strand","log2fc","input-count","output-count","log10pval","log10fdr","lfc-se","cigar","md-tag","project")
  printbed<-printbed[,c("chr","start","stop","ID","score","strand","log2FoldChange","ctrl_mean","exp_mean","pvalue","padj","lfcSE","project")]
  colnames(printbed)<-c("chr","start","stop","id","score","strand","log2fc","input-count","output-count","log10pval","log10fdr","lfc-se","project")
  printbed$strand[printbed$strand=="fwd"]="+"
  printbed$strand[printbed$strand=="rev"]="-"
  printbed$log10pval=-log10(printbed$log10pval)
  printbed$log10fdr=-log10(printbed$log10fdr)
  
  cat("printbed_In_function_1\n")
  cat(str(printbed))
  cat("\n")
  
  # write.table(printbed,paste0("results/",file_prefix,"_",'RNA',"_",fileDate(),".bed"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  write.table(printbed,paste0("results/",file_prefix,"_",'RNA',".bed"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  
  
  return(c(full_output, dds_results))
}


run_function = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  path_results<-paste(out,'results','/',sep='')
  
  if (file.exists(path_results)){
  }else{
    
    dir.create(file.path(path_results))
  }
  
  path_plots<-paste(out,'plots','/',sep='')
  
  if (file.exists(path_plots)){
  }else{
    
    dir.create(file.path(path_plots))
  }
  
  setwd(out)
  
  #### READ and transform Cell_Type_sel ----
  
  Tag_sel = opt$Tag_sel
  
  cat("Tag_sel_0\n")
  cat(sprintf(as.character(Tag_sel)))
  cat("\n")
  
  
  Cell_Type_sel<-gsub("__.+$","",Tag_sel)
  
  cat("Cell_Type_sel_0\n")
  cat(sprintf(as.character(Cell_Type_sel)))
  cat("\n")
  
  comparison_sel<-gsub("__.+$","",Tag_sel)
  
  cat("comparison_sel_0\n")
  cat(sprintf(as.character(comparison_sel)))
  cat("\n")
  
  #### Read countsData ----
  
  countsData <- readRDS(opt$countsData)
  
  
  cat("countsData_0\n")
  cat(str(countsData))
  cat("\n")
  
  #### Read attributesData ----
  
  attributesData <- readRDS(opt$attributesData)
  
  
  cat("attributesData_0\n")
  cat(str(attributesData))
  cat("\n")
  cat(str(unique(attributesData$Oligo)))
  cat("\n")
  
  #### Read attributesData_single_variants ----
  
  attributesData_single_variants <- readRDS(opt$attributesData_single_variants)
  
  
  cat("attributesData_single_variants_0\n")
  cat(str(attributesData_single_variants))
  cat("\n")
  cat(str(unique(attributesData_single_variants$Oligo)))
  cat("\n")
  
  attributesData<-rbind(attributesData,
                        attributesData_single_variants)
  
  cat("attributesData_1\n")
  cat(str(attributesData))
  cat("\n")
  cat(str(unique(attributesData$Oligo)))
  cat("\n")
  
  
  attributesData$ctrl_exp<-attributesData$Project
  
  indx.positive_controls<-grep("expCtrl",attributesData$Project)
  indx.negative_controls<-grep("negCtrl",attributesData$Project)
  
  attributesData$ctrl_exp[indx.positive_controls]<-'expCtrl'
  attributesData$ctrl_exp[indx.negative_controls]<-'negCtrl'
  
  cat("attributesData_1\n")
  cat(str(attributesData))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(attributesData$Project))))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(attributesData$ctrl_exp))))))
  cat("\n")
  
  colnames(attributesData)[which(colnames(attributesData) == 'Window')]<-'window'
  colnames(attributesData)[which(colnames(attributesData) == 'Allele')]<-'allele'
  colnames(attributesData)[which(colnames(attributesData) == 'Project')]<-'project'
  
  
  #### Read conditions_data_extended ----
  
  conditions_data_extended<-as.data.frame(fread(file=opt$conditions_data_extended, sep="\t", header=F), stringsAsFactors=F)
  
  cat("conditions_data_extended_0\n")
  cat(str(conditions_data_extended))
  cat("\n")
  
  
  #### Subset for Cell_Type_sel ----
  
  
  countsData_Oligo_sel<-countsData[which(countsData$Oligo%in%attributesData$ID),]
  
  
  cat("countsData_Oligo_sel_0\n")
  cat(str(countsData_Oligo_sel))
  cat("\n")
  
  if(dim(countsData_Oligo_sel)[1] >0){
    
    conditions_data<-conditions_data_extended[,-1]
    
    row.names(conditions_data)<-conditions_data_extended[,1]
    
    cat("conditions_data_0\n")
    cat(str(conditions_data))
    cat("\n")
    
    
    
    ##### test -----
    
    dataOut(countsData_Oligo_sel, conditions_data, attributesData, exclList = c(), altRef = T, Tag_sel, method = 'ss',negCtrlName="negCtrl",
            tTest=F, DEase=T, cSkew=T, correction="BH", cutoff=0.01, upDisp=T, prior=F, paired=F)
    
    
  }else{
    
    cat("No_oligos_in_attributesData_after_filtering\n")
    
  }#dim(countsData_Oligo_sel)[1] >0
  
  
 
}




printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--countsData"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--attributesData"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--attributesData_single_variants"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--conditions_data_extended"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Tag_sel"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  run_function(opt)
    
  
}


###########################################################################

system.time( main() )
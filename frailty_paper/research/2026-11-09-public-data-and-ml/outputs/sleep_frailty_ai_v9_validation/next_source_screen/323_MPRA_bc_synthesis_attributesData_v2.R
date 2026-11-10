
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
suppressMessages(library("Biostrings", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)


data_wrangling_design_equivalence = function(option_list)
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
  
 
  
  #### READ and transform multivariate_tracking ----
  
  multivariate_tracking = unlist(strsplit(opt$multivariate_tracking, split=";"))
  
  cat("multivariate_tracking\n")
  cat(sprintf(as.character(multivariate_tracking)))
  cat("\n")
  
  #### Read input equivalence_file ----
  
  equivalence_file<-readRDS(file = opt$equivalence_file)
  
  
  cat("equivalence_file_0\n")
  cat(str(equivalence_file))
  cat("\n")
 
  #### Read input design file ----
  
  design_file<-as.data.frame(fread(file = opt$design_file, sep="\t", header =F), stringsAsFactors=F)
  
  
  cat("design_file_0\n")
  cat(str(design_file))
  cat("\n")
  
  colnames(design_file)<-c("Error","Real_Tile","Carried_variants","proposed_3","seq_name","Label","factor4","Variant","chr","start","stop","VAR","PCR_5_arm","REST","REST_no_bc","bc","Enzymes","seq_name_bc","check")
  
  
  cat("design_file_1\n")
  cat(str(design_file))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(design_file$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(design_file$Label)))))
  cat("\n")
  
  ### Remove enzyme sequence and PCR 3' arms (rev comp of the oligo sequences)
  
  design_file$sequence<-gsub("^GGATCCGGTACC","",design_file$REST_no_bc)
  design_file$sequence<-gsub("CCTGCAGGGAATTC$|CACTGCGGCTCCTG$|GGTGCTCGCTATCG$","",design_file$sequence)
  # design_file$sequence<-gsub("CACTGCGGCTCCTG$","",design_file$sequence)
  # design_file$sequence<-gsub("GGTGCTCGCTATCG$","",design_file$sequence)
  
  
  design_file$sequence_length<-nchar(design_file$sequence)
  
  
  cat("design_file_2\n")
  cat(str(design_file))
  cat("\n")
  cat(sprintf(as.character(names(summary(design_file$sequence_length)))))
  cat("\n")
  cat(sprintf(as.character(summary(design_file$sequence_length))))
  cat("\n")
  
  design_file$KEY<-gsub(";.+","",design_file$seq_name_bc)
  
  
  cat("design_file_3\n")
  cat(str(design_file))
  cat("\n")
 
  
  design_file$original_TILE<-gsub("^[^;]+;","",design_file$seq_name_bc)
  design_file$original_TILE<-gsub(";.+","",design_file$original_TILE)
  
  design_file$TILE<-"NA"
  
  design_file$TILE[which(design_file$original_TILE == "NINE")]<-"TILE_1"
  design_file$TILE[which(design_file$original_TILE == "TWO_THIRDS")]<-"TILE_2"
  design_file$TILE[which(design_file$original_TILE == "HALF")]<-"TILE_3"
  design_file$TILE[which(design_file$original_TILE == "ONE_THIRD")]<-"TILE_4"
  design_file$TILE[which(design_file$original_TILE == "A_TENTH")]<-"TILE_5"
  
  design_file$carried_variants<-gsub("^[^;]+;[^;]+;[^;]+;[^;]+;","",design_file$seq_name_bc)
  design_file$carried_variants<-gsub(";.+","",design_file$carried_variants)
  design_file$carried_variants<-gsub("\\|","__chr",design_file$carried_variants)
  design_file$carried_variants<-paste("chr",design_file$carried_variants,sep='')
  
  design_file$Barcode<-gsub("^[^;]+;[^;]+;[^;]+;[^;]+;[^;]+;","",design_file$seq_name_bc)
  design_file$Barcode<-gsub(";.+","",design_file$Barcode)
  
  
  
  design_file$ALLELE_string<-gsub("^[^;]+;[^;]+;","",design_file$seq_name_bc)
  design_file$ALLELE_string<-gsub("_.+$","",design_file$ALLELE_string)
  design_file$ALLELE_string<-gsub("_.+$","",design_file$ALLELE_string)
  
  
  design_file$ALLELE_string<-gsub(";[0-9]+$","",design_file$ALLELE_string)
  
  indx.REF<-grep("REF", design_file$seq_name_bc)
  
  if(length(indx.REF) >0){
    
    design_file$carried_variants[indx.REF]<-"REF"
    design_file$carried_variants[-indx.REF]<-design_file$carried_variants[-indx.REF]
    
  }#length(indx.REF) >0
  
  
  design_file$Oligo<-paste(design_file$KEY,design_file$TILE,design_file$carried_variants, sep="__")
  
  cat("design_file_4\n")
  cat(str(design_file))
  cat("\n")
  
  NCGR_extract<-unique(design_file$KEY[which(design_file$Label == 'Negative_Control_Genomic_Regions')])
  
  cat("NCGR_extract\n")
  cat(sprintf(as.character(NCGR_extract)))
  cat("\n")
  
  ASE_CTRL_extract<-unique(design_file$KEY[which(design_file$Label == 'POSITIVE_CTRL')])
  
  cat("ASE_CTRL_extract\n")
  cat(sprintf(as.character(ASE_CTRL_extract)))
  cat("\n")
  
  Enhancer_CTRL_extract<-unique(design_file$KEY[which(design_file$Label == 'UNDETERMINED_CTRL')])
  
  cat("Enhancer_CTRL_extract\n")
  cat(sprintf(as.character(Enhancer_CTRL_extract)))
  cat("\n")
  
  Kousik_extract<-unique(design_file$KEY[which(design_file$Label == 'Kousik_variant')])
  
  cat("Kousik_extract\n")
  cat(sprintf(as.character(Kousik_extract)))
  cat("\n")
  
  #### READ and transform NCGR ----
  
  NCGR = unlist(strsplit(opt$NCGR, split=","))
  
  cat("NCGR_\n")
  cat(sprintf(as.character(NCGR)))
  cat("\n")
  
  #### READ and transform ASE_CTRL ----
  
  ASE_CTRL = unlist(strsplit(opt$ASE_CTRL, split=","))
  
  cat("ASE_CTRL_\n")
  cat(sprintf(as.character(ASE_CTRL)))
  cat("\n")
  
  #### READ and transform Enhancer_CTRL ----
  
  Enhancer_CTRL = unlist(strsplit(opt$Enhancer_CTRL, split=","))
  
  cat("Enhancer_CTRL_\n")
  cat(sprintf(as.character(Enhancer_CTRL)))
  cat("\n")
  
  #################--------------
  
  design_file_minus_NCGR<-unique(design_file[-which(design_file$Label == 'Negative_Control_Genomic_Regions'),])
  
  
  cat("design_file_minus_NCGR_0\n")
  cat(str(design_file_minus_NCGR))
  cat("\n")
  
  indx.int<-which(colnames(design_file_minus_NCGR)%in%c("seq_name_bc","sequence"))
  
  cat("indx.int_0\n")
  cat(str(indx.int))
  cat("\n")
  
  design_file_minus_NCGR_subset<-unique(design_file_minus_NCGR[,indx.int])
  
  cat("design_file_minus_NCGR_subset\n")
  cat(str(design_file_minus_NCGR_subset))
  cat("\n")
  
  
  
  design_file_NCGR<-unique(design_file[which(design_file$Label == 'Negative_Control_Genomic_Regions'),])
  
  
  cat("design_file_NCGR_0\n")
  cat(str(design_file_NCGR))
  cat("\n")
  
  
  
  indx.overlap<-which(design_file_minus_NCGR$seq_name_bc%in%equivalence_file$seq_name_bc)
  
  cat("indx.overlap_0\n")
  cat(str(indx.overlap))
  cat("\n")
  
  not_overlap_1<-design_file_minus_NCGR[-indx.overlap,]
  
  
  cat("not_overlap_1_0\n")
  cat(str(not_overlap_1))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(not_overlap_1$KEY))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(not_overlap_1$KEY)))))
  cat("\n")
  
  indx.overlap2<-which(equivalence_file$seq_name_bc%in%design_file_minus_NCGR$seq_name_bc)
  
  cat("indx.overlap2_0\n")
  cat(str(indx.overlap2))
  cat("\n")
  
  not_overlap_2<-equivalence_file[-indx.overlap2,]
  
  
  cat("not_overlap_2_0\n")
  cat(str(not_overlap_2))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(not_overlap_2$KEY))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(not_overlap_2$KEY)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(not_overlap_2$carried_variants))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(not_overlap_2$carried_variants)))))
  cat("\n")
  
  # 
  # indx.subset<-which(colnames(design_file_minus_NCGR)%in%c("seq_name_bc","sequence"))
  # 
  # 
  # design_file_minus_NCGR<-unique(design_file_minus_NCGR[,indx.subset])
  # 
  # cat("design_file_minus_NCGR_0\n")
  # cat(str(design_file_minus_NCGR))
  # cat("\n")
  
  indx.NCGR_sequence<-which(colnames(design_file_NCGR)%in%c("KEY","TILE","sequence"))
  
  cat("indx.NCGR_sequence_0\n")
  cat(str(indx.NCGR_sequence))
  cat("\n")
  
  design_file_NCGR_subset<-unique(design_file_NCGR[,indx.NCGR_sequence])
  
  cat("design_file_NCGR_subset\n")
  cat(str(design_file_NCGR_subset))
  cat("\n")
  
  not_overlap_2<-merge(not_overlap_2,
                       design_file_NCGR_subset,
                       by=c("KEY","TILE"))
  
  cat("not_overlap_2_1\n")
  cat(str(not_overlap_2))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(not_overlap_2$KEY))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(not_overlap_2$KEY)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(not_overlap_2$carried_variants))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(not_overlap_2$carried_variants)))))
  cat("\n")
  
  setwd(out)
  
  saveRDS(not_overlap_2, file="NCGR_for_port.rds")
  
  
  
  equivalence_file<-merge(equivalence_file,
                          design_file_minus_NCGR_subset,
                     by=c("seq_name_bc"))
  
 

  cat("equivalence_file_1\n")
  cat(str(equivalence_file))
  cat("\n")
  
  equivalence_file<-rbind(not_overlap_2,equivalence_file)
  
  
  
  cat("equivalence_file_2\n")
  cat(str(equivalence_file))
  cat("\n")
  
  
  check_multivariate<-equivalence_file[grep(multivariate_tracking,equivalence_file$Oligo),]
  
  cat("check_multivariate_0\n")
  cat(str(check_multivariate))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check_multivariate$carried_variants))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(check_multivariate$carried_variants)))))
  cat("\n")
  
  
  
  
  equivalence_file_subset<-unique(equivalence_file[,c(which(colnames(equivalence_file) == 'KEY'),which(colnames(equivalence_file) == 'TILE'),which(colnames(equivalence_file) == 'carried_variants'),
                                      which(colnames(equivalence_file) == 'Oligo'),which(colnames(equivalence_file) == 'sequence'))])
  
  cat("equivalence_file_subset_0\n")
  cat(str(equivalence_file_subset))
  cat("\n")
  cat(str(unique(equivalence_file_subset$Oligo)))
  cat("\n")
  
  
  equivalence_file_subset$Project<-'MPRA_bc_synthesis'
  
  equivalence_file_subset$Project[which(equivalence_file_subset$KEY%in%NCGR)]<-paste('MPRA_bc_synthesis','negCtrl', sep="_")
  equivalence_file_subset$Project[which(equivalence_file_subset$KEY%in%ASE_CTRL)]<-paste('MPRA_bc_synthesis','expCtrl', sep="_")
  equivalence_file_subset$Project[which(equivalence_file_subset$KEY%in%Enhancer_CTRL)]<-paste('MPRA_bc_synthesis','expCtrl', sep="_")
  
  
 
  equivalence_file_subset$Window<-NA
  
  
  equivalence_file_subset$Window[which(equivalence_file_subset$TILE == 'TILE_1')]<-'TILE_1'
  equivalence_file_subset$Window[which(equivalence_file_subset$TILE == 'TILE_2')]<-'TILE_2'
  equivalence_file_subset$Window[which(equivalence_file_subset$TILE == 'TILE_3')]<-'TILE_3'
  equivalence_file_subset$Window[which(equivalence_file_subset$TILE == 'TILE_4')]<-'TILE_4'
  equivalence_file_subset$Window[which(equivalence_file_subset$TILE == 'TILE_5')]<-'TILE_5'
  
  equivalence_file_subset$Window<-factor(equivalence_file_subset$Window,
                                 levels=c('TILE_1','TILE_2','TILE_3','TILE_4','TILE_5'),
                                 ordered=T)
  
  equivalence_file_subset$strand<-'+'
  
  cat("equivalence_file_subset_1\n")
  cat(str(equivalence_file_subset))
  cat("\n")
  cat(str(unique(equivalence_file_subset$Oligo)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(equivalence_file_subset$Project))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(equivalence_file_subset$Project)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(equivalence_file_subset$Allele))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(equivalence_file_subset$Allele)))))
  cat("\n")
  
  ###################################################################### SAVE
  
  setwd(out)
  
  saveRDS(equivalence_file_subset, file="ALL_sequences_unsorted.rds")
  
 
}

sort_diplotypes = function(option_list)
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
  
  
  
  #### READ and transform multivariate_tracking ----
  
  multivariate_tracking = unlist(strsplit(opt$multivariate_tracking, split=";"))
  
  cat("multivariate_tracking\n")
  cat(sprintf(as.character(multivariate_tracking)))
  cat("\n")
  
  setwd(out)
  
  ALL_sequences_unsorted<-readRDS(file="ALL_sequences_unsorted.rds")
  
  cat("ALL_sequences_unsorted_0\n")
  cat(str(ALL_sequences_unsorted))
  cat("\n")
  
  ##### Select diplotypes -----
  
  indx.diplotypes<-grep("__",ALL_sequences_unsorted$carried_variants)
  
  cat("indx.diplotypes_0\n")
  cat(str(indx.diplotypes))
  cat("\n")
  
  diplotypes<-ALL_sequences_unsorted[indx.diplotypes,]
  
  
  cat("diplotypes_0\n")
  cat(str(diplotypes))
  cat("\n")
  cat(str(unique(diplotypes$Oligo)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(diplotypes$carried_variants))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(diplotypes$carried_variants)))))
  cat("\n")
  
  
  #### LOOP through diplotypes ----
  
  
  diplo_KEYS<-unique(diplotypes$KEY)
  
  cat("diplo_KEYS_0\n")
  cat(str(diplo_KEYS))
  cat("\n")
  
  DEBUG<-0
  
  RESULT<-data.frame()
  for(i in 1:length(diplo_KEYS)){
    
    diplo_KEYS_sel<-diplo_KEYS[i]
    
    cat("------------------------------------>\t")
    cat(sprintf(as.character(diplo_KEYS_sel)))
    cat("\t")
    
    diplotypes_KEYS_sel<-diplotypes[which(diplotypes$KEY == diplo_KEYS_sel),]
    
    if(DEBUG ==1){
      
      cat("diplotypes_KEYS_sel_0\n")
      cat(str(diplotypes_KEYS_sel))
      cat("\n")
    }
    
    ALL_KEYS_sel<-ALL_sequences_unsorted[which(ALL_sequences_unsorted$KEY == diplo_KEYS_sel),]
    
    if(DEBUG ==1){
      
      cat("ALL_KEYS_sel_0\n")
      cat(str(ALL_KEYS_sel))
      cat("\n")
    }
    
    diplo_KEYS_Window<-unique(diplotypes_KEYS_sel$Window)
    
    for(k in 1:length(diplo_KEYS_Window)){
      
      diplo_KEYS_Window_sel<-diplo_KEYS_Window[k]
      
      # cat("------------------------------------>\n")
      cat(sprintf(as.character(diplo_KEYS_Window_sel)))
      cat("\t")
      
      
      diplotypes_KEYS_sel_Window_sel<-diplotypes_KEYS_sel[which(diplotypes_KEYS_sel$Window == diplo_KEYS_Window_sel),]
      
      if(DEBUG ==1){
        
        cat("diplotypes_KEYS_sel_Window_sel_0\n")
        cat(str(diplotypes_KEYS_sel_Window_sel))
        cat("\n")
      }
      
      
      ######### some tiles have multiple combinations of diplotypes
      
      array_diplotypes_KEYS_sel_Window_sel_carried_variants<-unique(diplotypes_KEYS_sel_Window_sel$carried_variants)
      
      
      if(DEBUG ==1){
        
        cat("array_diplotypes_KEYS_sel_Window_sel_carried_variants_0\n")
        cat(str(array_diplotypes_KEYS_sel_Window_sel_carried_variants))
        cat("\n")
      }
      
      for(m in 1:length(array_diplotypes_KEYS_sel_Window_sel_carried_variants)){
        
        array_diplotypes_KEYS_sel_Window_sel_carried_variants_sel<-array_diplotypes_KEYS_sel_Window_sel_carried_variants[m]
        
        # cat("------------------------------------>\n")
        cat(sprintf(as.character(array_diplotypes_KEYS_sel_Window_sel_carried_variants_sel)))
        cat("\t")
        
        
        diplotypes_KEYS_sel_Window_sel_carried_variants_sel<-diplotypes_KEYS_sel_Window_sel[which(diplotypes_KEYS_sel_Window_sel$carried_variants == array_diplotypes_KEYS_sel_Window_sel_carried_variants_sel),]
        
        if(DEBUG ==1){
          
          cat("diplotypes_KEYS_sel_Window_sel_carried_variants_sel_0\n")
          cat(str(diplotypes_KEYS_sel_Window_sel_carried_variants_sel))
          cat("\n")
        }
        
        ALL_KEYS_sel_Window_sel<-ALL_KEYS_sel[which(ALL_KEYS_sel$Window == diplo_KEYS_Window_sel),]
        
        if(DEBUG ==1){
          
          cat("ALL_KEYS_sel_Window_sel_0\n")
          cat(str(ALL_KEYS_sel_Window_sel))
          cat("\n")
        }
        
        REF_element<-ALL_KEYS_sel_Window_sel[which(ALL_KEYS_sel_Window_sel$carried_variants == "REF"),]
        
        if(DEBUG ==1){
          
          cat("REF_element_0\n")
          cat(str(REF_element))
          cat("\n")
          cat(sprintf(as.character(REF_element$carried_variants)))
          cat("\n")
        }
        
        
        
        variants_array<-unique(unlist(strsplit(diplotypes_KEYS_sel_Window_sel_carried_variants_sel$carried_variants, split= "__")))
        
        if(DEBUG ==1){
          
          cat("variants_array_0\n")
          cat(str(variants_array))
          cat("\n")
        }
        
        
        df<-data.frame()
        
        chr_vector<-NULL
        pos_vector<-NULL
        ref_vector<-NULL
        alt_vector<-NULL
        
        for(l in 1:length(variants_array)){
          
          variants_array_sel<-variants_array[l]
          
          # cat("------------------------------------>\n")
          cat(sprintf(as.character(variants_array_sel)))
          cat("\n")
          
          chr_vector[l]<-gsub("_.+$","",variants_array_sel)
          pos_vector[l]<-gsub("^[^_]+_","",variants_array_sel)
          pos_vector[l]<-as.integer(gsub("_.+$","",pos_vector[l]))
          ref_vector[l]<-gsub("^[^_]+_[^_]+_","",variants_array_sel)
          ref_vector[l]<-gsub("_.+$","",ref_vector[l])
          alt_vector[l]<-gsub("^[^_]+_[^_]+_[^_]+_","",variants_array_sel)
          
          tmp_element<-ALL_KEYS_sel_Window_sel[which(ALL_KEYS_sel_Window_sel$carried_variants == variants_array_sel),]
          
          allele_string<-rep("ref",length(variants_array))
          
          allele_string[l]<-'alt'
          
          tmp_element$Allele<-paste(allele_string, collapse='/')
          tmp_element$SNP<-diplotypes_KEYS_sel_Window_sel_carried_variants_sel$carried_variants
          
          
          if(DEBUG ==1){
            
            cat("tmp_element_0\n")
            cat(str(tmp_element))
            cat("\n")
            cat(str(chr_vector))
            cat("\n")
            cat(str(pos_vector))
            cat("\n")
            cat(str(ref_vector))
            cat("\n")
            cat(str(alt_vector))
            cat("\n")
          }
          
          
          
          df<-rbind(tmp_element,df)
          
          
        }#l in 1:length(variants_array)
        
        if(DEBUG ==1){
          
          cat("df_0\n")
          cat(str(df))
          cat("\n")
        }
        
        diplotypes_KEYS_sel_Window_sel_carried_variants_sel$Allele<-paste(rep('alt',length(variants_array)), collapse='/')
        diplotypes_KEYS_sel_Window_sel_carried_variants_sel$SNP<-diplotypes_KEYS_sel_Window_sel_carried_variants_sel$carried_variants
        
        
        
        REF_element$Allele<-paste(rep('ref',length(variants_array)), collapse='/')
        REF_element$SNP<-diplotypes_KEYS_sel_Window_sel_carried_variants_sel$carried_variants
        
        
        FINAL<-rbind(REF_element,df,diplotypes_KEYS_sel_Window_sel_carried_variants_sel)
        # FINAL$Project<-paste(FINAL$Project,"multiplotype",FINAL$SNP, sep="__")
        
        
        if(DEBUG ==1){
          
          cat("FINAL_0\n")
          cat(str(FINAL))
          cat("\n")
        }
        
        indx.REF_adaptation<-which(FINAL$carried_variants == 'REF') ### PROBLEM LINE
        
        if(DEBUG ==1){
          
          cat("indx.REF_adaptation_0\n")
          cat(str(indx.REF_adaptation))
          cat("\n")
        }
        
        FINAL$carried_variants[-indx.REF_adaptation]<-array_diplotypes_KEYS_sel_Window_sel_carried_variants_sel ### PROBLEM LINE
        
        if(DEBUG ==1){
          
          cat("FINAL_1\n")
          cat(str(FINAL))
          cat("\n")
        }
        
        FINAL$chr<-unique(chr_vector)
        FINAL$pos<-paste(pos_vector, collapse="__")
        FINAL$ref<-paste(ref_vector, collapse="__")
        FINAL$alt<-paste(alt_vector, collapse="__")
        
        if(DEBUG ==1){
          
          cat("FINAL_1\n")
          cat(str(FINAL))
          cat("\n")
        }
        
        RESULT<-rbind(FINAL,RESULT)
        
      }#m in 1:length(array_diplotypes_KEYS_sel_Window_sel_carried_variants)
    }#k in 1:length(diplo_KEYS_Window
  }#i in 1:length(diplo_KEYS)
  
  cat("RESULT_0\n")
  cat(str(RESULT))
  cat("\n")
  
  
  ###### SAVE RESULT --------------------
  
  unlink("Diplotypes.rds")
  
  setwd(out)
  
  saveRDS(RESULT, file="Diplotypes.rds")

  
}


sort_single_variants = function(option_list)
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
  
  #### READ and transform NCGR ----
  
  NCGR = unlist(strsplit(opt$NCGR, split=","))
  
  cat("NCGR_\n")
  cat(sprintf(as.character(NCGR)))
  cat("\n")
  
  #### READ and transform multivariate_tracking ----
  
  multivariate_tracking = unlist(strsplit(opt$multivariate_tracking, split=";"))
  
  cat("multivariate_tracking\n")
  cat(sprintf(as.character(multivariate_tracking)))
  cat("\n")
  
  setwd(out)
  
  ALL_sequences_unsorted<-readRDS(file="ALL_sequences_unsorted.rds")
  
  cat("ALL_sequences_unsorted_0\n")
  cat(str(ALL_sequences_unsorted))
  cat("\n")
  
  diplotypes<-readRDS(file="Diplotypes.rds")
  
  
  cat("diplotypes_0\n")
  cat(str(diplotypes))
  cat("\n")
  
  ##### Select single_variants -----
  
  diplotypes_carried_variants<-unique(diplotypes$carried_variants)
  
  cat("diplotypes_carried_variants_0\n")
  cat(str(diplotypes_carried_variants))
  cat("\n")
  
  diplotypes_carried_variants_minus_REF<-diplotypes_carried_variants[-which(diplotypes_carried_variants == 'REF')]
  
  cat("diplotypes_carried_variants_minus_REF_0\n")
  cat(str(diplotypes_carried_variants_minus_REF))
  cat("\n")
  
  diplotypes_carried_variants_minus_REF_broken_down_variants<-unique(unlist(strsplit(diplotypes_carried_variants_minus_REF, split="__")))
  
  cat("diplotypes_carried_variants_minus_REF_broken_down_variants_0\n")
  cat(str(diplotypes_carried_variants_minus_REF_broken_down_variants))
  cat("\n")
  
  
  single_variants<-ALL_sequences_unsorted[-which(ALL_sequences_unsorted$carried_variants%in%c(diplotypes_carried_variants_minus_REF,diplotypes_carried_variants_minus_REF_broken_down_variants)),]
  
  
  cat("single_variants_0\n")
  cat(str(single_variants))
  cat("\n")
  cat(str(unique(single_variants$Oligo)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(single_variants$carried_variants))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(single_variants$carried_variants)))))
  cat("\n")
  
  
  #### LOOP through single_variants ----
  
  
  single_KEYS<-unique(single_variants$KEY)
  
  cat("single_KEYS_0\n")
  cat(str(single_KEYS))
  cat("\n")
  
  DEBUG<-0
  
  RESULT<-data.frame()
  for(i in 1:length(single_KEYS)){
    
    single_KEYS_sel<-single_KEYS[i]
    
    cat("------------------------------------>\t")
    cat(sprintf(as.character(single_KEYS_sel)))
    cat("\t")
    
    FLAG_NCGR<-single_KEYS_sel[which(single_KEYS_sel%in%NCGR)]
    
    if(length(FLAG_NCGR) >0){
      
      DEBUG<-1
      
    }else{
      
      DEBUG<-0
      
    }#length(FLAG_NCGR) >0
    
    single_variants_KEYS_sel<-single_variants[which(single_variants$KEY == single_KEYS_sel),]
    
    
    
    
    if(DEBUG ==1){
      
      cat("single_variants_KEYS_sel_0\n")
      cat(str(single_variants_KEYS_sel))
      cat("\n")
    }
    
    
    single_KEYS_Window<-unique(single_variants_KEYS_sel$Window)
    
    for(k in 1:length(single_KEYS_Window)){
      
      single_KEYS_Window_sel<-single_KEYS_Window[k]
      
      # cat("------------------------------------>\n")
      cat(sprintf(as.character(single_KEYS_Window_sel)))
      cat("\t")
      
      
      single_variants_KEYS_sel_Window_sel<-single_variants_KEYS_sel[which(single_variants_KEYS_sel$Window == single_KEYS_Window_sel),]
      
      if(DEBUG ==1){
        
        cat("single_variants_KEYS_sel_Window_sel_0\n")
        cat(str(single_variants_KEYS_sel_Window_sel))
        cat("\n")
      }
      
      REF_element<-single_variants_KEYS_sel_Window_sel[which(single_variants_KEYS_sel_Window_sel$carried_variants == "REF"),]
      
      if(DEBUG ==1){
        
        cat("REF_element_0\n")
        cat(str(REF_element))
        cat("\n")
        cat(sprintf(as.character(REF_element$carried_variants)))
        cat("\n")
      }
      
      ######### some tiles have multiple combinations of diplotypes
      
      array_single_variants_KEYS_sel_Window_sel_carried_variants<-unique(single_variants_KEYS_sel_Window_sel$carried_variants)
      array_single_variants_KEYS_sel_Window_sel_carried_variants<-array_single_variants_KEYS_sel_Window_sel_carried_variants[which(array_single_variants_KEYS_sel_Window_sel_carried_variants != "REF")]
      
      if(DEBUG ==1){
        
        cat("array_single_variants_KEYS_sel_Window_sel_carried_variants_0\n")
        cat(str(array_single_variants_KEYS_sel_Window_sel_carried_variants))
        cat("\n")
      }
      
      if(length(array_single_variants_KEYS_sel_Window_sel_carried_variants) >0){
        
        for(m in 1:length(array_single_variants_KEYS_sel_Window_sel_carried_variants)){
          
          array_single_variants_KEYS_sel_Window_sel_carried_variants_sel<-array_single_variants_KEYS_sel_Window_sel_carried_variants[m]
          
          # cat("------------------------------------>\n")
          cat(sprintf(as.character(array_single_variants_KEYS_sel_Window_sel_carried_variants_sel)))
          cat("\t")
          
          
          single_variants_KEYS_sel_Window_sel_carried_variants_sel<-single_variants_KEYS_sel_Window_sel[which(single_variants_KEYS_sel_Window_sel$carried_variants == array_single_variants_KEYS_sel_Window_sel_carried_variants_sel),]
          
          if(DEBUG ==1){
            
            cat("single_variants_KEYS_sel_Window_sel_carried_variants_sel_0\n")
            cat(str(single_variants_KEYS_sel_Window_sel_carried_variants_sel))
            cat("\n")
          }
          
          variants_array<-unlist(strsplit(single_variants_KEYS_sel_Window_sel_carried_variants_sel$carried_variants, split= "__"))
          
          if(DEBUG ==1){
            
            cat("variants_array_0\n")
            cat(str(variants_array))
            cat("\n")
          }
          
          
          df<-data.frame()
          
          chr_vector<-NULL
          pos_vector<-NULL
          ref_vector<-NULL
          alt_vector<-NULL
          
          for(l in 1:length(variants_array)){
            
            variants_array_sel<-variants_array[l]
            
            # cat("------------------------------------>\n")
            cat(sprintf(as.character(variants_array_sel)))
            cat("\n")
            
            chr_vector[l]<-gsub("_.+$","",variants_array_sel)
            pos_vector[l]<-gsub("^[^_]+_","",variants_array_sel)
            pos_vector[l]<-as.integer(gsub("_.+$","",pos_vector[l]))
            ref_vector[l]<-gsub("^[^_]+_[^_]+_","",variants_array_sel)
            ref_vector[l]<-gsub("_.+$","",ref_vector[l])
            alt_vector[l]<-gsub("^[^_]+_[^_]+_[^_]+_","",variants_array_sel)
            
            tmp_element<-single_variants_KEYS_sel_Window_sel_carried_variants_sel[which(single_variants_KEYS_sel_Window_sel_carried_variants_sel$carried_variants == variants_array_sel),]
            
            allele_string<-rep("ref",length(variants_array))
            
            allele_string[l]<-'alt'
            
            tmp_element$Allele<-paste(allele_string, collapse='/')
            tmp_element$SNP<-single_variants_KEYS_sel_Window_sel_carried_variants_sel$carried_variants
            
            
            if(DEBUG ==1){
              
              cat("tmp_element_0\n")
              cat(str(tmp_element))
              cat("\n")
              cat(str(chr_vector))
              cat("\n")
              cat(str(pos_vector))
              cat("\n")
              cat(str(ref_vector))
              cat("\n")
              cat(str(alt_vector))
              cat("\n")
            }
            
            
            
            df<-rbind(tmp_element,df)
            
            
          }#l in 1:length(variants_array)
          
          if(DEBUG ==1){
            
            cat("df_0\n")
            cat(str(df))
            cat("\n")
          }
          
          REF_element$Allele<-paste(rep('ref',length(variants_array)), collapse='/')
          REF_element$SNP<-single_variants_KEYS_sel_Window_sel_carried_variants_sel$carried_variants
          
          
          FINAL<-rbind(REF_element,df)
          
          if(DEBUG ==1){
            
            cat("FINAL_0\n")
            cat(str(FINAL))
            cat("\n")
          }
          
          
          FINAL$chr<-unique(chr_vector)
          FINAL$pos<-paste(pos_vector, collapse="__")
          FINAL$ref<-paste(ref_vector, collapse="__")
          FINAL$alt<-paste(alt_vector, collapse="__")
          
          if(DEBUG ==1){
            
            cat("FINAL_1\n")
            cat(str(FINAL))
            cat("\n")
          }
          
          RESULT<-rbind(FINAL,RESULT)
          
          
          
        }#m in 1:length(array_single_variants_KEYS_sel_Window_sel_carried_variants)
        
      }else{
        
        # REF_element$Allele<-paste(rep('ref',length(variants_array)), collapse='/')
        # REF_element$SNP<-single_variants_KEYS_sel_Window_sel_carried_variants_sel$carried_variants
        # 
        # 
        # FINAL<-rbind(REF_element)
        # 
        # if(DEBUG ==1){
        #   
        #   cat("FINAL_0\n")
        #   cat(str(FINAL))
        #   cat("\n")
        # }
        # 
        # 
        # FINAL$chr<-unique(chr_vector)
        # FINAL$pos<-paste(pos_vector, collapse="__")
        # FINAL$ref<-paste(ref_vector, collapse="__")
        # FINAL$alt<-paste(alt_vector, collapse="__")
        # 
        # if(DEBUG ==1){
        #   
        #   cat("FINAL_1\n")
        #   cat(str(FINAL))
        #   cat("\n")
        # }
        # 
        # RESULT<-rbind(FINAL,RESULT)
        
        
        cat("----------------------------------------------------------------->\t")
        cat(sprintf(as.character("WARNING NO ALT ALLELE")))
        cat("\n")
        
        
      }#length(array_single_variants_KEYS_sel_Window_sel_carried_variants) >1
    }#k in 1:length(single_KEYS_Window

    # if(i ==2)
    # {
    #   ##########################################################################################################################################
    #   quit(status = 1)
    # }
  }#i in 1:length(single_KEYS)
  
  cat("RESULT_0\n")
  cat(str(RESULT))
  cat("\n")
  
  
  ###### SAVE RESULT --------------------
  
  
  setwd(out)
  
  unlink("Single_variants.rds")
  
  saveRDS(RESULT, file="Single_variants.rds")
  
}

annotationData_build_no_duplicates = function(option_list)
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
  
  #### READ and transform config_file ----
  
  config_file = as.data.frame(fread(file=opt$config_file, sep="\t", header = T), stringsAsFactors = F)
  
  cat("config_file_0\n")
  cat(str(config_file))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  path_sorted_attribute_data<-paste(out,'sorted_attribute_data','/',sep='')
  
  if (file.exists(path_sorted_attribute_data)){
    
  }else{
    dir.create(file.path(path_sorted_attribute_data))
  }
  
  ################### Read all the single varinats and the diplotypes -------------------
  
  setwd(out)

  single_variants<-readRDS(file="Single_variants.rds")
  
  
  cat("single_variants_0\n")
  cat(str(single_variants))
  cat("\n")
  
  ##### check duplicated Oligo in the single_variants
  
  indx.dup_single_variants<-which(duplicated(single_variants$Oligo) == "TRUE")
  
  cat("indx.dup_single_variants_0\n")
  cat(str(indx.dup_single_variants))
  cat("\n")
  
  
  if(length(indx.dup_single_variants) == 0){
    
    indx.int<-c(which(colnames(single_variants) == 'Oligo'),which(colnames(single_variants) == 'SNP'),which(colnames(single_variants) == 'Project'),
                which(colnames(single_variants) == 'Window'),which(colnames(single_variants) == 'strand'),which(colnames(single_variants) == 'Allele'),
                which(colnames(single_variants) == 'chr'),which(colnames(single_variants) == 'pos'),which(colnames(single_variants) == 'ref'),which(colnames(single_variants) == 'alt'))
    
    single_variants_subset<-unique(single_variants[,indx.int])
    
    colnames(single_variants_subset)[which(colnames(single_variants_subset) == 'Oligo')]<-'ID'
    
    cat("single_variants_subset_0\n")
    cat(str(single_variants_subset))
    cat("\n")
    cat(str(unique(single_variants_subset$ID)))
    cat("\n")
    
    setwd(path_sorted_attribute_data)
    
    saveRDS(single_variants_subset, file="attributesData_single_variants.rds")
    
  }else{
    
    stop("Duplicate Oligos in the single_variants file")
    
  }#length(indx.dup_single_variants) == 0
  
  setwd(out)
  
  diplotypes<-readRDS(file="Diplotypes.rds")
  
  
  cat("diplotypes_0\n")
  cat(str(diplotypes))
  cat("\n")
 
  
  
  #### MPRAmodel cannot take duplicate ids for Oligos
  #### However: 1. Diplotypes share the same ref allele comparison across the three diplotype groups
                # 2. Some diplotypes have multiple pairs of variants
  
  
  # first print the diplotype oligos in three groups based on the comparisons to do
  
  array_allele_string<-c('alt/ref','ref/alt','alt/alt')
  
  DEBUG<-1
  
  Files_att<-data.frame()
  
  for(i in 1:length(array_allele_string)){
    
    allele_sel<-array_allele_string[i]
    
    allele_comparison_name<-paste(gsub('/','_','ref/ref'), gsub('/','_',allele_sel),sep="_vs_")
      
      
    cat("--------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(allele_sel)))
    cat("\t")
    cat(sprintf(as.character(allele_comparison_name)))
    cat("\n")
    
    diplotypes_sel<-diplotypes[which(diplotypes$Allele%in%c('ref/ref',allele_sel)),]
    
    
    
    if(DEBUG == 1){
      
      cat("diplotypes_sel_0\n")
      cat(str(diplotypes_sel))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(diplotypes_sel$Allele))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(diplotypes_sel$Allele)))))
      cat("\n")
      
    }
    

    diplotypes_sel$Allele[which(diplotypes_sel$Allele == 'ref/ref')]<-'ref'
    diplotypes_sel$Allele[which(diplotypes_sel$Allele == allele_sel)]<-'alt'
    
    
    
    if(DEBUG == 1){
      
      cat("diplotypes_sel_1\n")
      cat(str(diplotypes_sel))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(diplotypes_sel$Allele))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(diplotypes_sel$Allele)))))
      cat("\n")
      
    }
   
    indx.dup<-which(duplicated(diplotypes_sel$Oligo) == "TRUE")
    
    if(DEBUG == 1){
      
      cat("indx.dup_0\n")
      cat(str(indx.dup))
      cat("\n")
    }
    
    if(length(indx.dup) >0){
      
      diplotypes_sel_NO_DUP<-diplotypes_sel[-indx.dup,]
      
      if(DEBUG == 1){
        
        cat("diplotypes_sel_NO_DUP_0\n")
        cat(str(diplotypes_sel_NO_DUP))
        cat("\n")
        cat(sprintf(as.character(names(summary(as.factor(diplotypes_sel_NO_DUP$Allele))))))
        cat("\n")
        cat(sprintf(as.character(summary(as.factor(diplotypes_sel_NO_DUP$Allele)))))
        cat("\n")
      }
      
      diplotypes_sel_DUP<-diplotypes_sel[indx.dup,]
      
      if(DEBUG == 1){
        
        cat("diplotypes_sel_DUP_0\n")
        cat(str(diplotypes_sel_DUP))
        cat("\n")
        cat(sprintf(as.character(names(summary(as.factor(diplotypes_sel_DUP$Allele))))))
        cat("\n")
        cat(sprintf(as.character(summary(as.factor(diplotypes_sel_DUP$Allele)))))
        cat("\n")
      }
      
      ## Second round
      
      indx.dup.2<-which(duplicated(diplotypes_sel_DUP$Oligo) == "TRUE")
      
      if(DEBUG == 1){
        
        cat("indx.dup.2_0\n")
        cat(str(indx.dup.2))
        cat("\n")
      }
      
      if(length(indx.dup.2) >0){
        
        diplotypes_sel_DUP_NO_DUP<-diplotypes_sel_DUP[-indx.dup.2,]
        
        if(DEBUG == 1){
          
          cat("diplotypes_sel_DUP_NO_DUP_0\n")
          cat(str(diplotypes_sel_DUP_NO_DUP))
          cat("\n")
          cat(sprintf(as.character(names(summary(as.factor(diplotypes_sel_DUP_NO_DUP$Allele))))))
          cat("\n")
          cat(sprintf(as.character(summary(as.factor(diplotypes_sel_DUP_NO_DUP$Allele)))))
          cat("\n")
        }
        
        diplotypes_sel_DUP_DUP<-diplotypes_sel_DUP[indx.dup.2,]
        
        if(DEBUG == 1){
          
          cat("diplotypes_sel_DUP_DUP_0\n")
          cat(str(diplotypes_sel_DUP_DUP))
          cat("\n")
          cat(sprintf(as.character(names(summary(as.factor(diplotypes_sel_DUP_DUP$Allele))))))
          cat("\n")
          cat(sprintf(as.character(summary(as.factor(diplotypes_sel_DUP_DUP$Allele)))))
          cat("\n")
        }
        
        ## Third round
        
        indx.dup.3<-which(duplicated(diplotypes_sel_DUP_DUP$Oligo) == "TRUE")
        
        if(DEBUG == 1){
          
          cat("indx.dup.3_0\n")
          cat(str(indx.dup.3))
          cat("\n")
        }
        
        #################################################### SAVE #########################################
        
        if(length(indx.dup.3) == 0){
          
          indx.int<-c(which(colnames(diplotypes_sel_NO_DUP) == 'Oligo'),which(colnames(diplotypes_sel_NO_DUP) == 'SNP'),which(colnames(diplotypes_sel_NO_DUP) == 'Project'),
                      which(colnames(diplotypes_sel_NO_DUP) == 'Window'),which(colnames(diplotypes_sel_NO_DUP) == 'strand'),which(colnames(diplotypes_sel_NO_DUP) == 'Allele'),
                      which(colnames(diplotypes_sel_NO_DUP) == 'chr'),which(colnames(diplotypes_sel_NO_DUP) == 'pos'),which(colnames(diplotypes_sel_NO_DUP) == 'ref'),which(colnames(diplotypes_sel_NO_DUP) == 'alt'))
          
          diplotypes_sel_NO_DUP_subset<-unique(diplotypes_sel_NO_DUP[,indx.int])
          
          colnames(diplotypes_sel_NO_DUP_subset)[which(colnames(diplotypes_sel_NO_DUP_subset) == 'Oligo')]<-'ID'
          
          if(DEBUG == 1){
            cat("diplotypes_sel_NO_DUP_subset_0\n")
            cat(str(diplotypes_sel_NO_DUP_subset))
            cat("\n")
            cat(str(unique(diplotypes_sel_NO_DUP_subset$ID)))
            cat("\n")
          }
          
          setwd(path_sorted_attribute_data)
          
          saveRDS(diplotypes_sel_NO_DUP_subset, file=paste("attributesData_diplotypes_1_",allele_comparison_name,".rds",sep=''))
          
          ##
          
          indx.int<-c(which(colnames(diplotypes_sel_DUP_NO_DUP) == 'Oligo'),which(colnames(diplotypes_sel_DUP_NO_DUP) == 'SNP'),which(colnames(diplotypes_sel_DUP_NO_DUP) == 'Project'),
                      which(colnames(diplotypes_sel_DUP_NO_DUP) == 'Window'),which(colnames(diplotypes_sel_DUP_NO_DUP) == 'strand'),which(colnames(diplotypes_sel_DUP_NO_DUP) == 'Allele'),
                      which(colnames(diplotypes_sel_DUP_NO_DUP) == 'chr'),which(colnames(diplotypes_sel_DUP_NO_DUP) == 'pos'),which(colnames(diplotypes_sel_DUP_NO_DUP) == 'ref'),which(colnames(diplotypes_sel_DUP_NO_DUP) == 'alt'))
          
          diplotypes_sel_DUP_NO_DUP_subset<-unique(diplotypes_sel_DUP_NO_DUP[,indx.int])
          
          colnames(diplotypes_sel_DUP_NO_DUP_subset)[which(colnames(diplotypes_sel_DUP_NO_DUP_subset) == 'Oligo')]<-'ID'
          
          if(DEBUG == 1){
            cat("diplotypes_sel_DUP_NO_DUP_subset_0\n")
            cat(str(diplotypes_sel_DUP_NO_DUP_subset))
            cat("\n")
            cat(str(unique(diplotypes_sel_DUP_NO_DUP_subset$ID)))
            cat("\n")
          }
          
          setwd(path_sorted_attribute_data)
          
          saveRDS(diplotypes_sel_DUP_NO_DUP_subset, file=paste("attributesData_diplotypes_2_",allele_comparison_name,".rds",sep=''))
          
          ##
          
          indx.int<-c(which(colnames(diplotypes_sel_DUP_DUP) == 'Oligo'),which(colnames(diplotypes_sel_DUP_DUP) == 'SNP'),which(colnames(diplotypes_sel_DUP_DUP) == 'Project'),
                      which(colnames(diplotypes_sel_DUP_DUP) == 'Window'),which(colnames(diplotypes_sel_DUP_DUP) == 'strand'),which(colnames(diplotypes_sel_DUP_DUP) == 'Allele'),
                      which(colnames(diplotypes_sel_DUP_DUP) == 'chr'),which(colnames(diplotypes_sel_DUP_DUP) == 'pos'),which(colnames(diplotypes_sel_DUP_DUP) == 'ref'),which(colnames(diplotypes_sel_DUP_DUP) == 'alt'))
          
          diplotypes_sel_DUP_DUP_subset<-unique(diplotypes_sel_DUP_DUP[,indx.int])
          
          colnames(diplotypes_sel_DUP_DUP_subset)[which(colnames(diplotypes_sel_DUP_DUP_subset) == 'Oligo')]<-'ID'
          
          if(DEBUG == 1){
            cat("diplotypes_sel_DUP_DUP_subset_0\n")
            cat(str(diplotypes_sel_DUP_DUP_subset))
            cat("\n")
            cat(str(unique(diplotypes_sel_DUP_DUP_subset$ID)))
            cat("\n")
          }
          
          setwd(path_sorted_attribute_data)
          
          saveRDS(diplotypes_sel_DUP_DUP_subset, file=paste("attributesData_diplotypes_3_",allele_comparison_name,".rds",sep=''))
          
          ##
          
          tmp_df<-as.data.frame(cbind(
            rbind(paste("diplotypes_1_",allele_comparison_name,sep=''),
                  paste("diplotypes_2_",allele_comparison_name,sep=''),
                  paste("diplotypes_3_",allele_comparison_name,sep='')),
            rbind(paste("attributesData_diplotypes_1_",allele_comparison_name,".rds",sep=''),
                  paste("attributesData_diplotypes_2_",allele_comparison_name,".rds",sep=''),
                  paste("attributesData_diplotypes_3_",allele_comparison_name,".rds",sep=''))), stringsAsFactors=F)
          
          colnames(tmp_df)<-c('comparison','attributesData')
          
          
          if(DEBUG == 1){
            cat("tmp_df_0\n")
            cat(str(tmp_df))
            cat("\n")
          }
          
          Files_att<-rbind(tmp_df,Files_att)
          
        }else{
          
          stop("Duplicate Oligos in the diplotypes file after 2 iterations")
          
        }#length(indx.dup.3) == 0
        
      }else{
        
        indx.int<-c(which(colnames(diplotypes_sel_NO_DUP) == 'Oligo'),which(colnames(diplotypes_sel_NO_DUP) == 'SNP'),which(colnames(diplotypes_sel_NO_DUP) == 'Project'),
                    which(colnames(diplotypes_sel_NO_DUP) == 'Window'),which(colnames(diplotypes_sel_NO_DUP) == 'strand'),which(colnames(diplotypes_sel_NO_DUP) == 'Allele'),
                    which(colnames(diplotypes_sel_NO_DUP) == 'chr'),which(colnames(diplotypes_sel_NO_DUP) == 'pos'),which(colnames(diplotypes_sel_NO_DUP) == 'ref'),which(colnames(diplotypes_sel_NO_DUP) == 'alt'))
        
        diplotypes_sel_NO_DUP_subset<-unique(diplotypes_sel_NO_DUP[,indx.int])
        
        colnames(diplotypes_sel_NO_DUP_subset)[which(colnames(diplotypes_sel_NO_DUP_subset) == 'Oligo')]<-'ID'
        
        if(DEBUG == 1){
          cat("diplotypes_sel_NO_DUP_subset_0\n")
          cat(str(diplotypes_sel_NO_DUP_subset))
          cat("\n")
          cat(str(unique(diplotypes_sel_NO_DUP_subset$ID)))
          cat("\n")
        }
        
        setwd(path_sorted_attribute_data)
        
        saveRDS(diplotypes_sel_NO_DUP_subset, file=paste("attributesData_diplotypes_1_",allele_comparison_name,".rds",sep=''))
        
        ###
        
        indx.int<-c(which(colnames(diplotypes_sel_DUP) == 'Oligo'),which(colnames(diplotypes_sel_DUP) == 'SNP'),which(colnames(diplotypes_sel_DUP) == 'Project'),
                    which(colnames(diplotypes_sel_DUP) == 'Window'),which(colnames(diplotypes_sel_DUP) == 'strand'),which(colnames(diplotypes_sel_DUP) == 'Allele'),
                    which(colnames(diplotypes_sel_DUP) == 'chr'),which(colnames(diplotypes_sel_DUP) == 'pos'),which(colnames(diplotypes_sel_DUP) == 'ref'),which(colnames(diplotypes_sel_DUP) == 'alt'))
        
        diplotypes_sel_DUP_subset<-unique(diplotypes_sel_DUP[,indx.int])
        
        colnames(diplotypes_sel_DUP_subset)[which(colnames(diplotypes_sel_DUP_subset) == 'Oligo')]<-'ID'
        
        if(DEBUG == 1){
          cat("diplotypes_sel_DUP_subset_0\n")
          cat(str(diplotypes_sel_DUP_subset))
          cat("\n")
          cat(str(unique(diplotypes_sel_DUP_subset$ID)))
          cat("\n")
        }
        
        setwd(path_sorted_attribute_data)
        
        saveRDS(diplotypes_sel_DUP_subset, file=paste("attributesData_diplotypes_2_",allele_comparison_name,".rds",sep=''))
        
        ##
        
        tmp_df<-as.data.frame(cbind(
          rbind(paste("diplotypes_1_",allele_comparison_name,sep=''),
                paste("diplotypes_2_",allele_comparison_name,sep='')),
          rbind(paste("attributesData_diplotypes_1_",allele_comparison_name,".rds",sep=''),
                paste("attributesData_diplotypes_2_",allele_comparison_name,".rds",sep=''))), stringsAsFactors=F)
        
        colnames(tmp_df)<-c('comparison','attributesData')
        
        
        if(DEBUG == 1){
          cat("tmp_df_0\n")
          cat(str(tmp_df))
          cat("\n")
        }
        
        Files_att<-rbind(tmp_df,Files_att)
        
      }#length(indx.dup.2) >0
      
      
    }else{
      
      indx.int<-c(which(colnames(diplotypes_sel) == 'Oligo'),which(colnames(diplotypes_sel) == 'SNP'),which(colnames(diplotypes_sel) == 'Project'),
                  which(colnames(diplotypes_sel) == 'Window'),which(colnames(diplotypes_sel) == 'strand'),which(colnames(diplotypes_sel) == 'Allele'),
                  which(colnames(diplotypes_sel) == 'chr'),which(colnames(diplotypes_sel) == 'pos'),which(colnames(diplotypes_sel) == 'ref'),which(colnames(diplotypes_sel) == 'alt'))
      
      diplotypes_sel_subset<-unique(diplotypes_sel[,indx.int])
      
      colnames(diplotypes_sel_subset)[which(colnames(diplotypes_sel_subset) == 'Oligo')]<-'ID'
      
      if(DEBUG == 1){
        cat("diplotypes_sel_subset_0\n")
        cat(str(diplotypes_sel_subset))
        cat("\n")
        cat(str(unique(diplotypes_sel_subset$ID)))
        cat("\n")
      }
      
      setwd(path_sorted_attribute_data)
      
      saveRDS(diplotypes_sel_subset, file=paste("attributesData_diplotypes_1_",allele_comparison_name,".rds",sep=''))
      
      ##
      
      tmp_df<-as.data.frame(cbind(
        rbind(paste("diplotypes_1_",allele_comparison_name,sep='')),
        rbind(paste("attributesData_diplotypes_1_",allele_comparison_name,".rds",sep=''))), stringsAsFactors=F)
      
      colnames(tmp_df)<-c('comparison','attributesData')
      
      
      if(DEBUG == 1){
        cat("tmp_df_0\n")
        cat(str(tmp_df))
        cat("\n")
      }
      
      Files_att<-rbind(tmp_df,Files_att)
      
      
    }#length(indx.dup) >0
    
   
    

    
  }#i in 1:length(array_allele_string
  
  cat("Files_att_0\n")
  cat(str(Files_att))
  cat("\n")
  
  # tmp_df<-as.data.frame(cbind(
  #   rbind(paste("single_variants",sep='')),
  #   rbind(paste("attributesData_single_variants.rds",sep=''))), stringsAsFactors=F)
  # 
  # colnames(tmp_df)<-c('comparison','attributesData')
  # 
  # Files_att<-rbind(tmp_df,Files_att)
  
  Files_att$attributesData<-paste(path_sorted_attribute_data,Files_att$attributesData, sep='')
  
  cat("Files_att_1\n")
  cat(str(Files_att))
  cat("\n")
  
  
  Final_config_file<-data.frame()
  
  DEBUG<-1
  
  for(i in 1:length(unique(config_file$Cell_Type))){
    
    Cell_Type_sel<-unique(config_file$Cell_Type)[i]
    
    cat("------------------------------------->\t")
    cat(sprintf(as.character(Cell_Type_sel)))
    cat("\n")
    
    tmp<-Files_att
    
    tmp$Cell_Type<-Cell_Type_sel
    
    if(DEBUG==1){
      
      cat("tmp_0\n")
      cat(str(tmp))
      cat("\n")
    }
    
    Final_config_file<-rbind(tmp,Final_config_file)
    
  }#i in 1:length(unique(config_file$Cell_Type)
  
  cat("Final_config_file_0\n")
  cat(str(Final_config_file))
  cat("\n")
  
  Final_config_file<-merge(config_file,
                           Final_config_file,
                           by="Cell_Type")
  
  cat("Final_config_file_1\n")
  cat(str(Final_config_file))
  cat("\n")
  
  setwd(out)
  
  write.table(Final_config_file, file="Run_file_countsData_and_condData_filtered_plus_attributesData.tsv", row.names=F, quote=F, sep="\t")
  
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
    make_option(c("--config_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--design_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--equivalence_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--multivariate_tracking"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--NCGR"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ASE_CTRL"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Enhancer_CTRL"), type="character", default=NULL, 
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

  data_wrangling_design_equivalence(opt)
  sort_diplotypes(opt)
  sort_single_variants(opt)
  annotationData_build_no_duplicates(opt)
 
  
}


###########################################################################

system.time( main() )
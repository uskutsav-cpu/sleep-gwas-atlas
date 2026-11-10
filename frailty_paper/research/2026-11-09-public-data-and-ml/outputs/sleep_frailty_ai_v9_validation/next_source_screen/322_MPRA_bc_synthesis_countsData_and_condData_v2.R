
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

oligoIsolate <- function(countsData, file_prefix){
  if("Error" %in% colnames(countsData)){
    countsData <- countsData[,c(1,2,7:dim(countsData)[2])]
  }
  tag_counts <- aggregate(. ~Oligo, data=countsData[,-1], FUN = sum)
  counts_oligo <- tag_counts[,-1]
  rownames(counts_oligo) <- tag_counts[,1]
  write.table(counts_oligo, paste0(file_prefix, "_counts.out"), quote = F, sep="\t")
  return(counts_oligo)
}

file_reader_to_countsData_LONG = function(option_list)
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
  
  #### READ and transform Kousik ----
  
  Kousik = unlist(strsplit(opt$Kousik, split=","))
  
  cat("Kousik_\n")
  cat(sprintf(as.character(Kousik)))
  cat("\n")
  
  #### READ and transform NCGR ----
  
  NCGR = unlist(strsplit(opt$NCGR, split=","))
  
  cat("NCGR_\n")
  cat(sprintf(as.character(NCGR)))
  cat("\n")
  
  
  #### READ and transform indir ----
  
  indir = opt$indir
  
  cat("indir_\n")
  cat(sprintf(as.character(indir)))
  cat("\n")
  
  setwd(indir)
  
  file_list <- list.files(path=indir, include.dirs = FALSE)
  
  cat("file_list\n")
  cat(str(file_list))
  cat("\n")
  
  indexes_sel <- grep("_gDNA_deduplicated_sorted_unique\\.tsv\\.gz$|_cDNA_deduplicated_sorted_unique\\.tsv\\.gz$",file_list)
  
  cat("indexes_sel\n")
  cat(sprintf(as.character(indexes_sel)))
  cat("\n")
  
  file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
  
  colnames(file_list_sel)<-"file"
  
  file_list_sel$file_extraction<-file_list_sel$file
  
  file_list_sel$file_extraction<-gsub("_gfpp","",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("ALL","K562",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("-GFP","R0minus",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("pGFP","R0plus",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("THP-1","THP1",file_list_sel$file_extraction)
  file_list_sel$file_extraction<-gsub("^R","K562_",file_list_sel$file_extraction)
  
  
  
  
  cat("file_list_sel_1\n")
  cat(str(file_list_sel))
  cat("\n")
  
  file_list_sel$Cell_Type<-gsub("_.+$","",file_list_sel$file_extraction)
  
  
  cat("file_list_sel_2\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
  file_list_sel$Replicate<-gsub("^[^_]+_","",file_list_sel$file_extraction)
  file_list_sel$Replicate<-gsub("_.+$","",file_list_sel$Replicate)
  
  #file_list_sel$Replicate[which(file_list_sel$Cell_Type == "ALL")]<-paste("ALL",file_list_sel$Replicate[which(file_list_sel$Cell_Type == "ALL")],sep='_')
  
  
  cat("file_list_sel_3\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
  file_list_sel$type<-gsub("^[^_]+_[^_]+_","",file_list_sel$file_extraction)
  file_list_sel$type<-gsub("_.+$","",file_list_sel$type)
  
  cat("file_list_sel_4\n")
  cat(str(file_list_sel))
  cat("\n")
  
  #file_list_sel$Cell_Type[which(file_list_sel$Cell_Type == "ALL")]<-"K562"
  
  file_list_sel$Sample_name<-paste(file_list_sel$Cell_Type,
                                   file_list_sel$Replicate,
                                   file_list_sel$type,
                              sep='_')
  
  file_list_sel<-file_list_sel[,-which(colnames(file_list_sel) == "file_extraction")]
  
  
  
  file_list_sel$type<-factor(file_list_sel$type,
                             c("gDNA","cDNA"),
                             ordered=T)
  file_list_sel$Cell_Type<-factor(file_list_sel$Cell_Type,
                                  c("K562","CHRF","HL60","THP1"),
                                  ordered=T)
  
  
  file_list_sel[order(file_list_sel$Cell_Type,file_list_sel$Replicate,file_list_sel$type),]
  
  
  
  
  cat("file_list_sel\n")
  cat(str(file_list_sel))
  cat("\n")
  cat(sprintf(as.character(names(summary(file_list_sel$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(file_list_sel$Cell_Type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(file_list_sel$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(file_list_sel$type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(file_list_sel$Replicate))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(file_list_sel$Replicate)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(file_list_sel$Sample_name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(file_list_sel$Sample_name)))))
  cat("\n")
  
 
  
  
  #### LOOP READ .gz files ----
  
  List_RESULTS<-list()
  
  List_equivalence<-list()
  
  # Results_DEF<-data.frame()
  
  DEBUG<-1
  
  for(i in 1:dim(file_list_sel)[1])
  {
    sel_file.df<-file_list_sel[i,]
    
    # cat("sel_file.df\n")
    # cat(str(sel_file.df))
    # cat("\n")
    
    setwd(indir)
    
    sel_file<-sel_file.df$file
    
    cat("------------------>\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(sel_file)))
    cat("\n")
    
    SIZE_gate<-file.info(sel_file)$size
    
    if(DEBUG ==1)
    {
      cat("SIZE_gate\n")
      cat(str(SIZE_gate))
      cat("\n")
    }
   
    
    
    if(SIZE_gate> 0)
    {
      
      
      LINE_gate<-length(readLines(sel_file))
      
      if(DEBUG ==1)
      {
        cat("LINE_gate\n")
        cat(str(LINE_gate))
        cat("\n")
      }
      
      if(LINE_gate> 0)
      {
        
        
        
        df<-as.data.frame(fread(file=sel_file, sep=" ", header = T, fill=TRUE), stringsAsFactors = F)
        colnames(df)<-c("counts","seq_name_bc")
        
        if(DEBUG ==1)
        {
          cat("df_0\n")
          cat(str(df))
          cat("\n")
        }
        
        df$Sample_name<-sel_file.df$Sample_name
        df$type<-sel_file.df$type
        df$Cell_Type<-sel_file.df$Cell_Type
        df$Replicate<-sel_file.df$Replicate
        
        df$KEY<-gsub(";.+","",df$seq_name_bc)
        
        if(DEBUG ==1)
        {
          cat("df_1\n")
          cat(str(df))
          cat("\n")
        }

        df$original_TILE<-gsub("^[^;]+;","",df$seq_name_bc)
        df$original_TILE<-gsub(";.+","",df$original_TILE)
        
        df$TILE<-"NA"
        
        df$TILE[which(df$original_TILE == "NINE")]<-"TILE_1"
        df$TILE[which(df$original_TILE == "TWO_THIRDS")]<-"TILE_2"
        df$TILE[which(df$original_TILE == "HALF")]<-"TILE_3"
        df$TILE[which(df$original_TILE == "ONE_THIRD")]<-"TILE_4"
        df$TILE[which(df$original_TILE == "A_TENTH")]<-"TILE_5"
        
        df$carried_variants<-gsub("^[^;]+;[^;]+;[^;]+;[^;]+;","",df$seq_name_bc)
        df$carried_variants<-gsub(";.+","",df$carried_variants)
        df$carried_variants<-gsub("\\|","__chr",df$carried_variants)
        df$carried_variants<-paste("chr",df$carried_variants,sep='')
        
        df$Barcode<-gsub("^[^;]+;[^;]+;[^;]+;[^;]+;[^;]+;","",df$seq_name_bc)
        df$Barcode<-gsub(";.+","",df$Barcode)
        
        
        
        df$ALLELE_string<-gsub("^[^;]+;[^;]+;","",df$seq_name_bc)
        df$ALLELE_string<-gsub("_.+$","",df$ALLELE_string)
        df$ALLELE_string<-gsub("_.+$","",df$ALLELE_string)
        

        df$ALLELE_string<-gsub(";[0-9]+$","",df$ALLELE_string)

        indx.REF<-grep("REF", df$seq_name_bc)
        
        if(length(indx.REF) >0){
          
          df$carried_variants[indx.REF]<-"REF"
          df$carried_variants[-indx.REF]<-df$carried_variants[-indx.REF]
          
        }#length(indx.REF) >0
        
        
        
        
        df$Oligo<-paste(df$KEY,df$TILE,df$carried_variants, sep="__")
        
        if(DEBUG ==1)
        {
          cat("df_2\n")
          cat(str(df))
          cat("\n")
        }
        
        
        indx.Kousik<-which(df$KEY%in%Kousik)
        
        if(DEBUG ==1)
        {
          cat("indx.Kousik_0\n")
          cat(str(indx.Kousik))
          cat("\n")
        }
        
        if(length(indx.Kousik) >0){
          
          df_minus_Kousik<-df[-indx.Kousik,]
          
        }else{
          
          df_minus_Kousik<-df
          
        }#length(indx.Kousik) >0
        
        if(DEBUG ==1)
        {
          cat("df_minus_Kousik_0\n")
          cat(str(df_minus_Kousik))
          cat("\n")
        }
        
          indx.NCGR<-which(df_minus_Kousik$KEY%in%NCGR)
          
          if(DEBUG ==1)
          {
            cat("indx.NCGR_0\n")
            cat(str(indx.NCGR))
            cat("\n")
          }
          
          if(length(indx.NCGR) >0){
            
            
            df_minus_Kousik_minus_NCGR<-df_minus_Kousik[-indx.NCGR,]
            
            if(DEBUG ==1)
            {
              cat("df_minus_Kousik_minus_NCGR_0\n")
              cat(str(df_minus_Kousik_minus_NCGR))
              cat("\n")
            }
            
            
            df_NCGR<-df[indx.NCGR,]
            
            df_NCGR$ALLELE_string<-gsub(paste(paste(rep(';NCGR', dim(df_NCGR)[1]),df_NCGR$Barcode, sep=';'), collapse = "|"),"",df_NCGR$ALLELE_string)
            
            
            if(DEBUG ==1)
            {
              cat("df_NCGR_0\n")
              cat(str(df_NCGR))
              cat("\n")
              cat(sprintf(as.character(names(summary(as.factor(df_NCGR$KEY))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(df_NCGR$KEY)))))
              cat("\n")
              
              cat(sprintf(as.character(names(summary(as.factor(df_NCGR$carried_variants))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(df_NCGR$carried_variants)))))
              cat("\n")
              
              cat(sprintf(as.character(names(summary(as.factor(df_NCGR$ALLELE_string))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(df_NCGR$ALLELE_string)))))
              cat("\n")
              
              cat(sprintf(as.character(names(summary(as.factor(df_NCGR$seq_name_bc))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(df_NCGR$seq_name_bc)))))
              cat("\n")
            }
            
            df_NCGR_prime<-df_NCGR
            
            
            for(iteration_NCGR in 1:length(NCGR)){
              
              NCGR_sel<-NCGR[iteration_NCGR]
              
              cat("------------------------------->\t")
              cat(sprintf(as.character(NCGR_sel)))
              cat("\t")
              
              mock_VAR<-NA
              
              if(NCGR_sel == 'Element_87'){
                
                mock_VAR<-'chr2_171570079_C_T'
                
              }else{
                if(NCGR_sel == 'Element_88'){
                  
                  mock_VAR<-'chr2_171580283_C_T'
                  
                }else{
                  
                  if(NCGR_sel == 'Element_89'){
                    
                    mock_VAR<-'chr2_175191238_G_A'
                    
                  }else{
                    
                    if(NCGR_sel == 'Element_90'){
                      
                      mock_VAR<-'chr2_175199760_T_C'
                      
                    }else{
                      
                    }#NCGR_sel == 'Element_90'
                    
                  }#NCGR_sel == 'Element_89'
                  
                }#NCGR_sel == 'Element_88'
                
              }#NCGR_sel == 'Element_87'
              
              cat(sprintf(as.character(mock_VAR)))
              cat("\t")
              
              indx.substitution<-which(df_NCGR_prime$KEY == NCGR_sel)
              
              if(DEBUG ==1)
              {
                cat("indx.substitution_0\n")
                cat(str(indx.substitution))
                cat("\n")
                
              }
              
              
              df_NCGR_prime$carried_variants[indx.substitution]<-mock_VAR
              df_NCGR_prime$ALLELE_string[indx.substitution]<-gsub("REF",mock_VAR,df_NCGR_prime$ALLELE_string[indx.substitution])
              df_NCGR_prime$seq_name_bc[indx.substitution]<-gsub("REF",mock_VAR,df_NCGR_prime$seq_name_bc[indx.substitution])
              
            }#iteration_NCGR in 1:length(NCGR)
            
            
            df_NCGR_prime$Oligo<-paste(df_NCGR_prime$KEY,df_NCGR_prime$TILE,df_NCGR_prime$carried_variants, sep="__")
            
            
            if(DEBUG ==1)
            {
              cat("df_NCGR_prime_0\n")
              cat(str(df_NCGR_prime))
              cat("\n")
              cat(sprintf(as.character(names(summary(as.factor(df_NCGR_prime$KEY))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(df_NCGR_prime$KEY)))))
              cat("\n")
              
              cat(sprintf(as.character(names(summary(as.factor(df_NCGR_prime$carried_variants))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(df_NCGR_prime$carried_variants)))))
              cat("\n")
              
              cat(sprintf(as.character(names(summary(as.factor(df_NCGR_prime$ALLELE_string))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(df_NCGR_prime$ALLELE_string)))))
              cat("\n")
              
              cat(sprintf(as.character(names(summary(as.factor(df_NCGR_prime$seq_name_bc))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(df_NCGR_prime$seq_name_bc)))))
              cat("\n")
            }
            
            
            
            df_minus_Kousik_minus_NCGR_plus_NCGR_prime<-rbind(df_minus_Kousik_minus_NCGR,df_NCGR,df_NCGR_prime)
            
            if(DEBUG ==1)
            {
              cat("df_minus_Kousik_minus_NCGR_REMEMBER\n")
              cat(str(df_minus_Kousik_minus_NCGR))
              cat("\n")
              
              cat("df_minus_Kousik_minus_NCGR_plus_NCGR_prime_0\n")
              cat(str(df_minus_Kousik_minus_NCGR_plus_NCGR_prime))
              cat("\n")
            
            }
            
            
            

          }else{
            
            df_minus_Kousik_minus_NCGR_plus_NCGR_prime<-df_minus_Kousik_minus_NCGR
            
            if(DEBUG ==1)
            {
              cat("df_minus_Kousik_minus_NCGR_REMEMBER\n")
              cat(str(df_minus_Kousik_minus_NCGR))
              cat("\n")
              
              cat("df_minus_Kousik_minus_NCGR_plus_NCGR_prime_WRONG\n")
              cat(str(df_minus_Kousik_minus_NCGR_plus_NCGR_prime))
              cat("\n")
              
            }
            
            

          }#length(indx.NCGR) >0
          
          
         
          
         
        
        
        df_subset_equivalence<-unique(df_minus_Kousik_minus_NCGR_plus_NCGR_prime[,c(which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'KEY'),which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'original_TILE'),which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'TILE'),
                                which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'ALLELE_string'), which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'carried_variants'), which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'Oligo'),which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'Barcode'),
                                which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'seq_name_bc'))])
        
        if(DEBUG ==1)
        {
          cat("df_subset_equivalence\n")
          cat(str(df_subset_equivalence))
          cat("\n")
        }
        
        List_equivalence[[i]]<-df_subset_equivalence
        
        
        df_subset<-unique(df_minus_Kousik_minus_NCGR_plus_NCGR_prime[,c(which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'Oligo'),which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'Barcode'),
                                which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'Sample_name'), which(colnames(df_minus_Kousik_minus_NCGR_plus_NCGR_prime) == 'counts'))])
        
        if(DEBUG ==1)
        {
          cat("df_subset\n")
          cat(str(df_subset))
          cat("\n")
        }
        
        List_RESULTS[[i]]<-df_subset
        
        # Results_DEF<-rbind(df_subset,Results_DEF)
        
        
      }#LINE_gate
      else{
        
        cat(sprintf("empty_file\n"))
        cat(sprintf(as.character(sel_file)))
        cat("\n")
        
        
      }
    }#SIZE_gates
  }# i files
  
 
  countsData_LONG = unique(as.data.frame(data.table::rbindlist(List_RESULTS, fill = T)))
  
  cat("countsData_LONG_0\n")
  cat(str(countsData_LONG))
  cat("\n")
  
  df_equivalence = unique(as.data.frame(data.table::rbindlist(List_equivalence, fill = T)))
  
  cat("df_equivalence_0\n")
  cat(str(df_equivalence))
  cat("\n")
  
  
  #### SAVE FILES ----
  
  filename=paste('countsData_LONG','.rds',sep='')
  
  setwd(out)
  saveRDS(countsData_LONG,file=filename)
  
  #### SAVE FILES ----
  
  filename=paste('df_equivalence','.rds',sep='')
  
  setwd(out)
  saveRDS(df_equivalence,file=filename)
}

filter_function = function(option_list)
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
  
  
  path_filtered_count_data<-paste(out,'filtered_count_data','/',sep='')
  
  if (file.exists(path_filtered_count_data)){
    
  }else{
    dir.create(file.path(path_filtered_count_data))
  }
  
  
  setwd(out)
  
  countsData_LONG<-readRDS(file="countsData_LONG.rds")
  
  cat("countsData_LONG_0\n")
  cat(str(countsData_LONG))
  cat("\n")
  cat(str(unique(countsData_LONG$Sample_name)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(countsData_LONG$Sample_name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(countsData_LONG$Sample_name)))))
  cat("\n")
  
  
  #### READ and transform filter_counts_gDNA ----
  
  filter_counts_gDNA = opt$filter_counts_gDNA
  
  cat("filter_counts_gDNA_\n")
  cat(sprintf(as.character(filter_counts_gDNA)))
  cat("\n")
  
  #### READ and transform filter_counts_cDNA ----
  
  filter_counts_cDNA = opt$filter_counts_cDNA
  
  cat("filter_counts_cDNA_\n")
  cat(sprintf(as.character(filter_counts_cDNA)))
  cat("\n")
  
  #### READ and transform Pseudocount ----
  
  Pseudocount = opt$Pseudocount
  
  cat("filter_counts_cDNA_\n")
  cat(sprintf(as.character(Pseudocount)))
  cat("\n")
  
  #### READ and transform Replicates_selected ----
  
  Replicates_selected = unlist(strsplit(opt$Replicates_selected, split=","))
  
  cat("Replicates_selected_\n")
  cat(str(Replicates_selected))
  cat("\n")
  
  Replicates_selected_gDNA = paste(Replicates_selected,"_gDNA",sep='')
  
  cat("Replicates_selected_gDNA_0\n")
  cat(str(Replicates_selected_gDNA))
  cat("\n")
  
  Replicates_selected_cDNA = paste(Replicates_selected,"_cDNA",sep='')
  
  cat("Replicates_selected_cDNA_0\n")
  cat(str(Replicates_selected_cDNA))
  cat("\n")
  
  #### Replicates accepted countsData_LONG ------------------
  
  countsData_LONG_replicates_accepted<-countsData_LONG[which(countsData_LONG$Sample_name%in%c(Replicates_selected_gDNA,Replicates_selected_cDNA)),]
  
  
  cat("countsData_LONG_replicates_accepted_0\n")
  cat(str(countsData_LONG_replicates_accepted))
  cat("\n")
  cat(str(unique(countsData_LONG_replicates_accepted$Sample_name)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(countsData_LONG_replicates_accepted$Sample_name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(countsData_LONG_replicates_accepted$Sample_name)))))
  cat("\n")
  

  countsData_LONG_replicates_accepted$type<-gsub("^[^_]+_[^_]+_","",countsData_LONG_replicates_accepted$Sample_name)
  countsData_LONG_replicates_accepted$Cell_Type<-gsub("_.+$","",countsData_LONG_replicates_accepted$Sample_name)
  countsData_LONG_replicates_accepted$Replicate<-gsub("^[^_]+_","",countsData_LONG_replicates_accepted$Sample_name)
  countsData_LONG_replicates_accepted$Replicate<-gsub("_.+$","",countsData_LONG_replicates_accepted$Replicate)
  
  cat("countsData_LONG_replicates_accepted_1\n")
  cat(str(countsData_LONG_replicates_accepted))
  cat("\n")
  cat(str(unique(countsData_LONG_replicates_accepted$Sample_name)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(countsData_LONG_replicates_accepted$type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(countsData_LONG_replicates_accepted$type)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(countsData_LONG_replicates_accepted$Cell_Type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(countsData_LONG_replicates_accepted$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(countsData_LONG_replicates_accepted$Replicate))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(countsData_LONG_replicates_accepted$Replicate)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(countsData_LONG_replicates_accepted$Sample_name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(countsData_LONG_replicates_accepted$Sample_name)))))
  cat("\n")
  
  countsData_LONG_replicates_accepted$type<-factor(countsData_LONG_replicates_accepted$type,
                             c("gDNA","cDNA"),
                             ordered=T)
  countsData_LONG_replicates_accepted$Cell_Type<-factor(countsData_LONG_replicates_accepted$Cell_Type,
                                  c("CHRF","K562","HL60","THP1"),
                                  ordered=T)
  
  

  cat("countsData_LONG_replicates_accepted_2\n")
  cat(str(countsData_LONG_replicates_accepted))
  cat("\n")
  cat(str(unique(countsData_LONG_replicates_accepted$Oligo)))
  cat("\n")
  cat(str(unique(countsData_LONG_replicates_accepted$Sample_name)))
  cat("\n")
  cat(sprintf(as.character(names(summary(countsData_LONG_replicates_accepted$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(countsData_LONG_replicates_accepted$type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(countsData_LONG_replicates_accepted$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(countsData_LONG_replicates_accepted$Cell_Type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(countsData_LONG_replicates_accepted$Replicate)))))
  cat("\n")
  cat(sprintf(as.character(summary(countsData_LONG_replicates_accepted$Replicate))))
  cat("\n")
  
  countsData_LONG_replicates_accepted.dt<-data.table(countsData_LONG_replicates_accepted, key=c("Cell_Type","Replicate","type","Sample_name","Oligo"))
  
  
  countsData_LONG_replicates_accepted_collapse_barcodes<-as.data.frame(countsData_LONG_replicates_accepted.dt[,.(Filter_counts=sum(counts)),by=key(countsData_LONG_replicates_accepted.dt)], stringsAsFactors=F)
  
  cat("countsData_LONG_replicates_accepted_collapse_barcodes_0\n")
  cat(str(countsData_LONG_replicates_accepted_collapse_barcodes))
  cat("\n")
  cat(str(unique(countsData_LONG_replicates_accepted_collapse_barcodes$Oligo)))
  cat("\n")
  cat(str(unique(countsData_LONG_replicates_accepted_collapse_barcodes$Sample_name)))
  cat("\n")
  cat(sprintf(as.character(names(summary(countsData_LONG_replicates_accepted_collapse_barcodes$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(countsData_LONG_replicates_accepted_collapse_barcodes$type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(countsData_LONG_replicates_accepted_collapse_barcodes$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(summary(countsData_LONG_replicates_accepted_collapse_barcodes$Cell_Type))))
  cat("\n")
  cat(sprintf(as.character(names(summary(countsData_LONG_replicates_accepted_collapse_barcodes$Replicate)))))
  cat("\n")
  cat(sprintf(as.character(summary(countsData_LONG_replicates_accepted_collapse_barcodes$Replicate))))
  cat("\n")
  
  ##### LOOP Cell types and filtered oligo by cell type -----
  
  DEBUG<-1
  
  run_master_file_df<-data.frame()
  
  for(i in 1:length(levels(countsData_LONG_replicates_accepted_collapse_barcodes$Cell_Type))){
    
    Cell_Type_sel<-levels(countsData_LONG_replicates_accepted_collapse_barcodes$Cell_Type)[i]
    
    cat("--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------->\n")
    cat("--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(Cell_Type_sel)))
    cat("\n")
    
    df_Cell_Type_sel<-droplevels(countsData_LONG_replicates_accepted_collapse_barcodes[which(countsData_LONG_replicates_accepted_collapse_barcodes$Cell_Type == Cell_Type_sel),])
    
    if(DEBUG == 1){
      cat("df_Cell_Type_sel_0\n")
      cat(str(df_Cell_Type_sel))
      cat("\n")
    }
    
    df_Cell_Type_sel$FILTER<-NA
    
    df_Cell_Type_sel$FILTER[which(df_Cell_Type_sel$type == 'gDNA' & df_Cell_Type_sel$Filter_counts >= filter_counts_gDNA)]<-'PASS'
    df_Cell_Type_sel$FILTER[which(df_Cell_Type_sel$type == 'gDNA' & df_Cell_Type_sel$Filter_counts < filter_counts_gDNA)]<-'FAIL'
    df_Cell_Type_sel$FILTER[which(df_Cell_Type_sel$type == 'cDNA' & df_Cell_Type_sel$Filter_counts >= filter_counts_cDNA)]<-'PASS'
    df_Cell_Type_sel$FILTER[which(df_Cell_Type_sel$type == 'cDNA' & df_Cell_Type_sel$Filter_counts < filter_counts_cDNA)]<-'FAIL'
    
    df_Cell_Type_sel$FILTER<-factor(df_Cell_Type_sel$FILTER,
                                    levels=c('FAIL','PASS'),
                                    ordered=T)
    
    if(DEBUG == 1){
      cat("df_Cell_Type_sel_1\n")
      cat(str(df_Cell_Type_sel))
      cat("\n")
      cat(sprintf(as.character(names(summary(df_Cell_Type_sel$FILTER)))))
      cat("\n")
      cat(sprintf(as.character(summary(df_Cell_Type_sel$FILTER))))
      cat("\n")
    }
    
    
    df_Cell_Type_sel.dt<-data.table(df_Cell_Type_sel, key=c("Cell_Type","Replicate","Oligo"))
    
    
    df_Cell_Type_sel_collapse_filters<-as.data.frame(df_Cell_Type_sel.dt[,.(FILTER=paste(unique(FILTER), collapse = ";")),by=key(df_Cell_Type_sel.dt)], stringsAsFactors=F)
    df_Cell_Type_sel_collapse_filters$FILTER<-factor(df_Cell_Type_sel_collapse_filters$FILTER)
    
    if(DEBUG == 1){
      cat("df_Cell_Type_sel_collapse_filters_0\n")
      cat(str(df_Cell_Type_sel_collapse_filters))
      cat("\n")
      cat(str(unique(df_Cell_Type_sel_collapse_filters$Oligo)))
      cat("\n")
      cat(sprintf(as.character(names(summary(df_Cell_Type_sel_collapse_filters$FILTER)))))
      cat("\n")
      cat(sprintf(as.character(summary(df_Cell_Type_sel_collapse_filters$FILTER))))
      cat("\n")
    }
    
    df_Cell_Type_sel_collapse_filters_PASS<-df_Cell_Type_sel_collapse_filters[which(df_Cell_Type_sel_collapse_filters$FILTER == 'PASS'),]
    
    if(DEBUG == 1){
      cat("df_Cell_Type_sel_collapse_filters_PASS_0\n")
      cat(str(df_Cell_Type_sel_collapse_filters_PASS))
      cat("\n")
      cat(str(unique(df_Cell_Type_sel_collapse_filters_PASS$Oligo)))
      cat("\n")
      cat(sprintf(as.character(names(summary(df_Cell_Type_sel_collapse_filters_PASS$FILTER)))))
      cat("\n")
      cat(sprintf(as.character(summary(df_Cell_Type_sel_collapse_filters_PASS$FILTER))))
      cat("\n")
    }
    
    
   Oligos_accepted_Cell_Type_sel<-unique(df_Cell_Type_sel_collapse_filters_PASS$Oligo)
   
   if(DEBUG == 1){
     cat("Oligos_accepted_Cell_Type_sel_0\n")
     cat(str(Oligos_accepted_Cell_Type_sel))
     cat("\n")
   }
   
   countsData_LONG_replicates_accepted_Cell_Type_sel<-droplevels(countsData_LONG_replicates_accepted[which(countsData_LONG_replicates_accepted$Cell_Type == Cell_Type_sel),])
   
   if(DEBUG == 1){
     cat("countsData_LONG_replicates_accepted_Cell_Type_sel_0\n")
     cat(str(countsData_LONG_replicates_accepted_Cell_Type_sel))
     cat("\n")
     cat(str(unique(countsData_LONG_replicates_accepted_Cell_Type_sel$Cell_Type)))
     cat("\n")
     cat(str(unique(countsData_LONG_replicates_accepted_Cell_Type_sel$Oligo)))
     cat("\n")
   }
   
   countsData_LONG_replicates_accepted_Cell_Type_sel_PASS<-countsData_LONG_replicates_accepted_Cell_Type_sel[which(countsData_LONG_replicates_accepted_Cell_Type_sel$Oligo%in%Oligos_accepted_Cell_Type_sel),]
   
   if(DEBUG == 1){
     cat("countsData_LONG_replicates_accepted_Cell_Type_sel_PASS_0\n")
     cat(str(countsData_LONG_replicates_accepted_Cell_Type_sel_PASS))
     cat("\n")
     cat(str(unique(countsData_LONG_replicates_accepted_Cell_Type_sel_PASS$Cell_Type)))
     cat("\n")
     cat(str(unique(countsData_LONG_replicates_accepted_Cell_Type_sel_PASS$Oligo)))
     cat("\n")
     cat(sprintf(as.character(names(summary(countsData_LONG_replicates_accepted_Cell_Type_sel_PASS$counts)))))
     cat("\n")
     cat(sprintf(as.character(summary(countsData_LONG_replicates_accepted_Cell_Type_sel_PASS$counts))))
     cat("\n")
   }
   
   countsData_LONG_replicates_accepted_Cell_Type_sel_PASS$counts<-countsData_LONG_replicates_accepted_Cell_Type_sel_PASS$counts + Pseudocount
   
   if(DEBUG == 1){
     cat("countsData_LONG_replicates_accepted_Cell_Type_sel_PASS_1_Pseudocount_added\n")
     cat(str(countsData_LONG_replicates_accepted_Cell_Type_sel_PASS))
     cat("\n")
     cat(str(unique(countsData_LONG_replicates_accepted_Cell_Type_sel_PASS$Cell_Type)))
     cat("\n")
     cat(str(unique(countsData_LONG_replicates_accepted_Cell_Type_sel_PASS$Oligo)))
     cat("\n")
     cat(sprintf(as.character(names(summary(countsData_LONG_replicates_accepted_Cell_Type_sel_PASS$counts)))))
     cat("\n")
     cat(sprintf(as.character(summary(countsData_LONG_replicates_accepted_Cell_Type_sel_PASS$counts))))
     cat("\n")
   }
   

   countsData_wide<-data.frame(pivot_wider(countsData_LONG_replicates_accepted_Cell_Type_sel_PASS,
                                           id_cols=c("Barcode","Oligo"),
                                           names_from=Sample_name,
                                           values_from=counts), stringsAsFactors=F)
   if(DEBUG == 1){
     cat("countsData_wide_0\n")
     cat(str(countsData_wide))
     cat("\n")
     cat(str(unique(countsData_wide$Oligo)))
     cat("\n")
   }
   
   countsData_wide[is.na(countsData_wide)]<-Pseudocount
   
   
   if(DEBUG == 1){
     cat("countsData_wide_1\n")
     cat(str(countsData_wide))
     cat("\n")
   }
   
   #### test oligoIsolate ----
   

   setwd(path_results)
   
   oligoIsolate(countsData_wide, Cell_Type_sel)
   
   #### create also conditions data ----
   
   
   Sample_names<-colnames(countsData_wide)[-which(colnames(countsData_wide)%in%c('Barcode','Oligo'))]
   
   if(DEBUG == 1){
     cat("Sample_names_0\n")
     cat(str(Sample_names))
     cat("\n")
   }
   
   #### create conditions data file -----
   
   
   conditions_data <- data.frame(matrix(vector(), 
                                        length(Sample_names), 
                                        2,
                                        dimnames=list(c(),
                                                      c("type","cell_type"))),
                                 stringsAsFactors=F)
   
   if(DEBUG == 1){
     cat("conditions_data_0\n")
     cat(str(conditions_data))
     cat("\n")
   }
   
   row.names(conditions_data)<-Sample_names
   
   conditions_data$type<-gsub("^[^_]+_[^_]+_","",row.names(conditions_data))
   conditions_data$cell_type<-gsub("_.+$","",row.names(conditions_data))
   
   conditions_data$type[which(conditions_data$type == "gDNA")]<-'DNA'
   conditions_data$type[which(conditions_data$type == "cDNA")]<-'RNA'
   
   
   if(DEBUG == 1){
     cat("conditions_data_1\n")
     cat(str(conditions_data))
     cat("\n")
   }
   
   #### save both count and condition data ----
   
   filename_counts=paste(Cell_Type_sel,'_countsData_wide_filtered','.rds',sep='')
   filename_cond=paste(Cell_Type_sel,'_condData_filtered','.tsv',sep='')
   
     
   setwd(path_filtered_count_data)
   
   saveRDS(countsData_wide,file=filename_counts)
   
   write.table(conditions_data, file=filename_cond, row.names=T, col.names = F, quote=F, sep="\t")
   
   temp_df<-as.data.frame(cbind(Cell_Type_sel,filename_counts,filename_cond), stringsAsFactors=F)
   colnames(temp_df)<-c("Cell_Type","countsData","condData")
   
   temp_df$countsData<-paste(path_filtered_count_data,temp_df$countsData,sep="")
   temp_df$condData<-paste(path_filtered_count_data,temp_df$condData,sep="")
   
   if(DEBUG == 1){
     cat("temp_df_0\n")
     cat(str(temp_df))
     cat("\n")
   }
   
   run_master_file_df<-rbind(temp_df,run_master_file_df)
    
  }#i in 1:length(levels(countsData_LONG_replicates_accepted_collapse_barcodes$Cell_type)
  
  
  cat("run_master_file_df_0\n")
  cat(str(run_master_file_df))
  cat("\n")
  
  
  #### SAVE MASTER FILE TO RUN MPRAmodel----
  
  setwd(out)
  
  write.table(run_master_file_df, file="Run_file_countsData_and_condData_filtered.tsv", row.names=F, quote=F, sep="\t")
  
  
  
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
    make_option(c("--indir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Replicates_selected"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--filter_counts_gDNA"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--filter_counts_cDNA"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--NCGR"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Kousik"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Pseudocount"), type="numeric", default=NULL, 
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
  
  file_reader_to_countsData_LONG(opt)
  filter_function(opt)
  
    
  
}


###########################################################################

system.time( main() )
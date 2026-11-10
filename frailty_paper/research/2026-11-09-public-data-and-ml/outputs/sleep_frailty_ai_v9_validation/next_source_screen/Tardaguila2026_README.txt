## Documentation for the manuscript: "Integrating Natural and Engineered Genetic Variation to Decode Regulatory Influence on Blood Traits"
Tardaguila et al, 2025

## 00_conda_environments_and_R_libraries/

   -Folder that contains the yml files to recreate the conda environments and to install the R libraries with the versions used in the paper

## 01_Tables/

   - Folder with the main and supplementary Tables of the paper.

     - Table_S2/ includes Jupyter_notebooks/ and Dependencies/ folders to recreate it
     - Table_S3/ includes Jupyter_notebooks/ and Dependencies/ folders to recreate it
     - Table_S4/ includes Jupyter_notebooks/ and Dependencies/ folders to recreate it
     - Table_S5/ includes Jupyter_notebooks/ and Dependencies/ folders to recreate it
     - Table_S6/ includes Jupyter_notebooks/ and Dependencies/ folders to recreate it


## 02_Figures/

   - Folder with the main and supplementary figures of the paper.

   All figure folders include Jupyter_notebooks/ to recreate each figure panel

## 03_data/

   - Folder with intermediate data points that allow the re-analysis by other users:

     -MPRA/NEW_pipeline_ALL_TOGETHER_analysis/ -> The gDNA and cDNA count matrices for all the MPRA replicates (*deduplicated_sorted_unique.tsv.gz)
     -GenIE/ -> the original rgenie plots for each of the variants assayed
     -Flow_Cyt_ILR_analysis/FlowCyt_results_for_model_added_all_quadrants_v2.csv -> proportions of each of the cell types in the Flow Cytometry used in the ILR model
     -rs139141690_CUX1/merged_clusters_final_annotated.rds -> Seurat object with all the cells that passed the QC and were genotyped with all the annotations. It is the starting point of the DE, DA, chromVAR and GRN analyses.

## 04_data_analysis/

   - Folder with all the scripts employed in the manuscript

     -Variant_scores/ -> Scripts to perform the statistical comparison across the prioritization tiers and with the known pathogenic variants
     -ukbiobank_fbc/ -> Scripts used for the conditional analysis on common variants
     -MPRA/Alignment_and_counts/ -> Scripts used to align the gDNA and cDNA MPRA fastq files and count the barcodes
     -MPRA/Analysis/ -> Scripts used to pre-process the counts (03_data/MPRA/NEW_pipeline_ALL_TOGETHER_analysis/) to the MPRA model format, run MPRA model, perform the meta-analysis across tiles for every variant in every cell type and finally the lasso regression to find predicitive variables of log2FC
     -Bulk_RNAseq/INTERVAL_analysis/ -> Scripts to performe the DE and ATU analysis on the INTERVAL whole blood dataset
     -Bulk_RNAseq/BluePrint_analysis/ -> Scripts to performe the DE and ATU analysis on the BluePrint immune cell isolates
     -Bulk_RNAseq/Active_Pathways_ORA_with_Dorothea/ -> Scripts to perform the ORA on the bulk RNA-seq DE results
     -Bulk_RNAseq/Bash_scripts/ -> Scripts to perform the Chisq test on ATU variants for the MPRA positive/ negative label, the test of the number of heterozygote carriers in the RNA no effect samples and finally the lasso regression for predictive variables of the MPRA-DE/ATU+ group
     -Bulk_RNAseq/rs35592432_FOXP1/ -> Scripts to perform DE and ATU (DRIMSeq model) in the rs35592432-FOXP1 clone lines
     -INTERVAL_sQTLs/ -> Scripts to find proxy variants to sQTL variants from INTERVAL in our group of selected genes
     -GTEx_sQTLs/ -> Scripts to find proxy variants to sQTL variants from GTEx in our group of selected genes
     -locus_plots/ -> Scripts to generate the locus plots including GWAS, DE/ATU, ATAC-seq, PCHi-C and MPRA results illustrated in Figure 4 and Supplementary Figure 4.
     -GenIE/ -> Scripts used in the analysis of the GenIE fastq libraries
     -TF_motif_prediction/All_motifs_in_MPRA/ -> Scripts used to predict TF motifs using JASPAR and HOMER for all the variants in Table S6
     -TF_motif_prediction/rs139141690_Del_80bp_and_Del_16bp/ -> Scripts used predict TF motifs using JASPAR and HOMER for rs139141690, the 16 bp and the 80 bp deletions
     -Flow_Cyt_ILR_analysis/ -> Code for the ILR modeling of the Flow Cytometry proportions (03_data/-Flow_Cyt_ILR_analysis/FlowCyt_results_for_model_added_all_quadrants_v2.csv)
     -multiome/rs139141690_K562/ -> Scripts used in the QC, genotyping of the multiome experiment. After QC and genotyping (03_data/rs139141690_CUX1/merged_clusters_final_annotated.rds) we include the scripts to perform the DE (plus ORA), DA, chromVAR and the GRN with SIMBA.

## 05_Add_hoc_analysis_for_referees/

   ATU_5_prime/ -> Exploration of TSS and variants for ATU cases
   Reanalysis_MPV_UKB_GWAS/ -> Re-analysis of MPV in UK-B WGS to search for rs139141690 hozygous alternative individuals
   Figures/ -> Figures for referee's questions

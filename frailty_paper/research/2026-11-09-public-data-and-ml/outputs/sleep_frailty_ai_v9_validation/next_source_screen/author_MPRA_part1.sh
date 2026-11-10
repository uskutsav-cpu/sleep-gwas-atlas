#!/bin/bash

eval "$(conda shell.bash hook)"
 
 
#### Save in a new folder, Analysis_filtered
   
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

Threshold_log2FC_meta=$3 # Changed. 0.255 .In https://github.com/julirsch/finemapped_mpra/blob/main/code/preprocess/mpra_meta.R is 1.
Threshold_FC_meta_padj=$4 # 0.01  Same as https://github.com/julirsch/finemapped_mpra/blob/main/code/preprocess/mpra_meta.R

Threshold_Skew_meta_padj=$(echo '0.1') # 0.1 Same as https://github.com/julirsch/finemapped_mpra/blob/main/code/preprocess/mpra_meta.R


echo "Threshold_log2FC_meta: $Threshold_log2FC_meta"" ""Threshold_FC_meta_padj: $Threshold_FC_meta_padj"



output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files

# ### build_counts_data


type=$(echo "build_counts_data")
outfile_build_counts_data=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_build_counts_data
echo -n "" > $outfile_build_counts_data
name_build_counts_data=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")


Rscript_build_counts_data=$(echo "$Rscripts_path""322_MPRA_bc_synthesis_countsData_and_condData_v2.R")

#indir=$(echo "/group/soranzo/manuel.tardaguila/MPRA_bc_synthesis/MPRA_bc_synthesis_counts/NEW_pipeline_ALL_TOGETHER_analysis/")
indir=$(echo "/group/soranzo/manuel.tardaguila/Zenodo_V2F_paper/03_data/MPRA/NEW_pipeline_ALL_TOGETHER_analysis/")
filter_counts_gDNA=$(echo '5')
filter_counts_cDNA=$(echo '5')
Pseudocount=$(echo '1')
Replicates_selected=$(echo 'K562_6,K562_7,K562_14,K562_16,K562_17,K562_18,K562_19,CHRF_R1,CHRF_R2,CHRF_R11,CHRF_R12,CHRF_R13,CHRF_R14,CHRF_R15,HL60_R5,HL60_R7,HL60_R9,HL60_R10,HL60_R12,HL60_R14,HL60_R15,THP1_R0plus,THP1_R1,THP1_R6,THP1_R7,THP1_R8,THP1_R9,THP1_R10')
Kousik=$(echo 'Kousik_1,Kousik_10,Kousik_11,Kousik_2,Kousik_3,Kousik_4,Kousik_5,Kousik_6,Kousik_7,Kousik_8,Kousik_9')
NCGR=$(echo 'Element_87,Element_88,Element_89,Element_90')

myjobid_build_counts_data=$(sbatch --output=$outfile_build_counts_data --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=8 --mem-per-cpu=4096 --parsable --job-name $name_build_counts_data --wrap="Rscript $Rscript_build_counts_data --indir $indir --filter_counts_gDNA $filter_counts_gDNA --filter_counts_cDNA $filter_counts_cDNA --Pseudocount $Pseudocount --Replicates_selected $Replicates_selected --Kousik $Kousik --NCGR $NCGR --type $type --out $output_dir")
myjobid_seff_build_counts_data=$(sbatch --dependency=afterany:$myjobid_build_counts_data --open-mode=append --output=$outfile_build_counts_data --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_build_counts_data >> $outfile_build_counts_data")

############################################################################################################################## build_attribute_data


type=$(echo "build_attribute_data")
outfile_build_attribute_data=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_build_attribute_data
echo -n "" > $outfile_build_attribute_data
name_build_attribute_data=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")
 

Rscript_build_attribute_data=$(echo "$Rscripts_path""323_MPRA_bc_synthesis_attributesData_v2.R")

design_file=$(echo "/group/soranzo/manuel.tardaguila/MPRA_bc_synthesis/Correspondence_barcode_variant.txt")
equivalence_file=$(echo "$output_dir""df_equivalence.rds")
multivariate_tracking=$(echo "Element_115__TILE_3")
NCGR=$(echo 'Element_87,Element_88,Element_89,Element_90')
ASE_CTRL=$(echo 'Element_1,Element_15,Element_16,Element_26,Element_42,Element_49,Element_64')
Enhancer_CTRL=$(echo 'Element_14,Element_21,Element_3,Element_30,Element_35,Element_51,Element_56,Element_57,Element_63,Element_65')
config_file=$(echo "$output_dir""Run_file_countsData_and_condData_filtered.tsv")

# --dependency=afterany:$myjobid_build_counts_data

myjobid_build_attribute_data=$(sbatch --dependency=afterany:$myjobid_build_counts_data --output=$outfile_build_attribute_data --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=1024 --parsable --job-name $name_build_attribute_data --wrap="Rscript $Rscript_build_attribute_data --design_file $design_file --equivalence_file $equivalence_file --multivariate_tracking $multivariate_tracking --NCGR $NCGR --ASE_CTRL $ASE_CTRL --Enhancer_CTRL $Enhancer_CTRL --config_file $config_file --type $type --out $output_dir")
myjobid_seff_build_attribute_data=$(sbatch --dependency=afterany:$myjobid_build_attribute_data --open-mode=append --output=$outfile_build_attribute_data --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_build_attribute_data >> $outfile_build_attribute_data")

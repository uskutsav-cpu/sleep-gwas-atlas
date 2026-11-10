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

config_file=$(echo "$output_dir""Run_file_countsData_and_condData_filtered_plus_attributesData.tsv")

echo "$config_file"

cd "$output_dir"

declare -a arr


end_of_file=0
counter=0
array=()
array_Merge=()

while [[ $end_of_file == 0 ]]
do
  read -r line
  end_of_file=$?


  ((counter++))

    echo "counter:$counter  LINE:$line"

  if [ $counter != 1 ]; then

      if [[ "$line" =~ [^[:space:]] ]]; then

        a=($(echo "$line" | tr "\t" '\n'))

        Cell_Type=${a[0]}
        comparison=${a[3]}
        countsData=${a[1]}
        condData=${a[2]}
	attributesData=${a[4]}


#        echo "Cell_Type:$Cell_Type  comparison:$comparison  countsData:$countsData condData:$condData attributesData:$attributesData"

	tag=$(echo "$Cell_Type""__""$comparison")

	echo "------------------------------------------------------------------------------------->$tag"

        # ####################################################################################################################################################################### run_MPRAmodel
	# ####################################################################################################################################################################### run_MPRAmodel
	# ####################################################################################################################################################################### run_MPRAmodel
	

	type=$(echo "$tag""_run_MPRAmodel")
	outfile_run_MPRAmodel=$(echo "$Log_files""outfile_3_""$type"".log")
	touch $outfile_run_MPRAmodel
	echo -n "" > $outfile_run_MPRAmodel
	name_run_MPRAmodel=$(echo "$type""_job")
	seff_name=$(echo "seff""_""$type")


	Rscript_run_MPRAmodel=$(echo "$Rscripts_path""324_MPRA_bc_synthesis_MPRA_model_adapted_v3.R")


	conditions_data_extended=$(echo "$condData")
	attributesData_single_variants=$(echo "$output_dir""sorted_attribute_data""/""attributesData_single_variants.rds")


	myjobid_run_MPRAmodel=$(sbatch --output=$outfile_run_MPRAmodel --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096 --parsable --job-name $name_run_MPRAmodel --wrap="Rscript $Rscript_run_MPRAmodel --countsData $countsData --attributesData $attributesData --conditions_data_extended $conditions_data_extended --Tag_sel $tag --attributesData_single_variants $attributesData_single_variants --type $type --out $output_dir")
	myjobid_seff_run_MPRAmodel=$(sbatch --dependency=afterany:$myjobid_run_MPRAmodel --open-mode=append --output=$outfile_run_MPRAmodel --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_run_MPRAmodel >> $outfile_run_MPRAmodel
")


	arr[${#arr[@]}]="$myjobid_run_MPRAmodel"



      fi #no spaces
    fi #counter > 1
done < "$config_file"

done_string=$(echo "--dependency=afterany:"""""${arr[@]}"""")
echo "$done_string"

dependency_string=$(echo $done_string|sed -r 's/ /:/g')

echo "$dependency_string"

####### ############################ collapse_and_meta_analysis
####### ############################ collapse_and_meta_analysis
####### ############################ collapse_and_meta_analysis
 

type=$(echo "collapse_and_meta_analysis")
outfile_collapse_and_meta_analysis=$(echo "$Log_files""outfile_4_""$type"".log")
touch $outfile_collapse_and_meta_analysis
echo -n "" > $outfile_collapse_and_meta_analysis
name_collapse_and_meta_analysis=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")


Rscript_collapse_and_meta_analysis=$(echo "$Rscripts_path""325_MPRA_bc_synthesis_Global_meta_analysis_v5.R")

NCGR=$(echo 'Element_87,Element_88,Element_89,Element_90')
ASE_CTRL=$(echo 'Element_1,Element_15,Element_16,Element_26,Element_42,Element_49,Element_64')
Enhancer_CTRL=$(echo 'Element_14,Element_21,Element_3,Element_30,Element_35,Element_51,Element_56,Element_57,Element_63,Element_65')


# $dependency_string

myjobid_collapse_and_meta_analysis=$(sbatch $dependency_string --output=$outfile_collapse_and_meta_analysis --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=1024 --parsable --job-name $name_collapse_and_meta_analysis --wrap="Rscript $Rscript_collapse_and_meta_analysis --Threshold_log2FC_meta $Threshold_log2FC_meta --Threshold_FC_meta_padj $Threshold_FC_meta_padj --Threshold_Skew_meta_padj $Threshold_Skew_meta_padj --NCGR $NCGR --ASE_CTRL $ASE_CTRL --Enhancer_CTRL $Enhancer_CTRL --type $type --out $output_dir")
myjobid_seff_collapse_and_meta_analysis=$(sbatch --dependency=afterany:$myjobid_collapse_and_meta_analysis --open-mode=append --output=$outfile_collapse_and_meta_analysis --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collapse_and_meta_analysis >> $outfile_collapse_and_meta_analysis")


#### Add_VAR #############################


type=$(echo "Add_VAR""_""$analysis")
outfile_Add_VAR=$(echo "$Log_files""outfile_5_""$type"".log")
touch $outfile_Add_VAR
echo -n "" > $outfile_Add_VAR
name_Add_VAR=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_Add_VAR=$(echo "$Rscripts_path""368_MPRA_endpoint_files_v4.R")

Table_S6=$(echo "/group/soranzo/manuel.tardaguila/Paper_bits/FIX_TABLES/Provisional_Tables/Table_S6_Provisional.rds")
MPRA_alignment_results="/group/soranzo/manuel.tardaguila/MPRA_explore_results/Alingment_check/MPRA_bc_synthesis_check_reconstruct_from_alignment.rds"
MPRA_result_SNP_and_CT="$output_dir""/collapsed_results/MPRA_results_meta_analysis_collapsed_by_SNP_and_Cell_Type_Threshold_log2FC_""$Threshold_log2FC_meta""_Threshold_FC_meta_padj_""$Threshold_FC_meta_padj"".rds"
MPRA_prior_to_meta_analysis="$output_dir""/collapsed_results/MPRA_results_prior_to_meta_analysis.tsv"

myjobid_Add_VAR=$(sbatch --dependency=afterany:$myjobid_collapse_and_meta_analysis --job-name=$name_Add_VAR --output=$outfile_Add_VAR --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="Rscript $Rscript_Add_VAR --Table_S6 $Table_S6 --MPRA_alignment_results $MPRA_alignment_results --MPRA_result_SNP_and_CT $MPRA_result_SNP_and_CT --MPRA_prior_to_meta_analysis $MPRA_prior_to_meta_analysis --type $type --out $output_dir")
myjobid_seff_Add_VAR=$(sbatch --dependency=afterany:$myjobid_Add_VAR --open-mode=append --output=$outfile_Add_VAR --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Add_VAR >> $outfile_Add_VAR")



ll ## 1. Get the regular matrix, add by barcode and normalize

$ bash Bash_scripts/95_MPRA_bc_synthesis_MPRAmodel_v3_PartONE.sh /group/soranzo/manuel.tardaguila/Zenodo_V2F_paper/04_data_analysis/MPRA/Analysis/ Analysis_filtered_Plus_Pseudocount 0.255 0.01

## 2. Perform MPRA model and meta-analysis

$ bash Bash_scripts/95_MPRA_bc_synthesis_MPRAmodel_v3_PartTWO.sh /group/soranzo/manuel.tardaguila/Zenodo_V2F_paper/04_data_analysis/MPRA/Analysis/ Analysis_filtered_Plus_Pseudocount 0.255 0.01

## 3. Perform lasso regression

$ bash Bash_scripts/114_paper_relaunch_MPRA_v4_lasso_only_cell_matched.sh /group/soranzo/manuel.tardaguila/Zenodo_V2F_paper/04_data_analysis/MPRA/Analysis/Analysis_filtered_Plus_Pseudocount/ Lasso_regression

  ## Data from K562 RNA from 2024-12-02 11:37:05.333693      https://www.encodeproject.org/rnaget-report/?type=RNAExpression&file.assay_title=polyA+plus+RNA-seq&file.biosample_ontology.classification=cell+line&file.biosample_ontology.term_name=K562&limit=all

## 3. All plots

$ bash ~/Scripts/Wraper_scripts/115_paper_relaunch_MPRA_initial_plots_v2.sh /group/soranzo/manuel.tardaguila/paper_relaunch/MPRA_S_no_global_metanalysis/ Initial_plots

## 3. Overlap with previous studies

$ bash ~/Scripts/Wraper_scripts/112_Explore_overlap_with_Ulirsch.sh /group/soranzo/manuel.tardaguila/MPRA_explore_results/ Explore_overlap_with_Ulirsch /group/soranzo/manuel.tardaguila/MPRA_bc_synthesis/Analysis_filtered_Plus_Pseudocount/collapsed_results/MPRA_results_meta_analysis_collapsed_by_SNP_Threshold_log2FC_0.255_Threshold_FC_meta_padj_0.01.rds

## 4. Figure panels:

-------------> Jupyter notebooks: Figure_2_panel_B.ipynb
-------------> Jupyter notebooks: Figure_2_panel_C.ipynb
-------------> Jupyter notebooks: Figure_2_panel_D.ipynb
-------------> Jupyter notebooks: Figure_2_panel_E.ipynb
-------------> Jupyter notebooks: Figure_S1_panel_C.ipynb
-------------> Jupyter notebooks: Figure_S1_panel_D.ipynb
-------------> Jupyter notebooks: Figure_S1_panel_E.ipynb
-------------> Jupyter notebooks: Figure_S1_panel_F.ipynb
-------------> Jupyter notebooks: Figure_S1_panel_G.ipynb




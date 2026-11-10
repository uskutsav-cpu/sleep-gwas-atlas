# MPRA_bc_synthesis_Sample_alignment_and_counts Sample Preparation


############# SAMPLE PROCESSING


command line:
 
$ bash Sample_processing.sh sample_file_GLOBAL.txt \<path to leave output files\> \<batch\> \<path to the MPRA_programmed_pipeline\> \<path to the original fastq files produced by bcl2fastq\> \<mem\> \<processors\> \<bsub queue\>

example:

$ bash Sample_processing.sh examples/sample_file_example_gDNA.txt Output_folder/ 141 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/141/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/ 20000 5 normal

$ nohup bash ~/Scripts/MPRA_programmed_pipeline/Sample_processing.sh /nfs/team151_data03/MPRA_Programmed/sample_file_GLOBAL.txt /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/HT_TEST/ HT_TEST /nfs/users/nfs_m/mt19/Scripts/MPRA_programmed_pipeline/ /nfs/team151_data03/MPRA_Programmed/ 16000 4 normal &

$ bash Sample_processing.sh examples/sample_file_GLOBAL.txt Output_folder/ ALL_TOGETHER /nfs/team151_data03/MPRA_Programmed/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/MPRA_Programmed_pipeline/ 20000 5 normal


Memory requirements: 16-48 GB  (for R16_gDNA)

##  carries the information of how to merge the fastqs and where are the UMIs and the R1 and R2

$ head examples/sample_file_example.txt

CHRF_R1_gDNA_S4_L001_R1_001.fastq.gz    CHRF_R1_gDNA_S4_L002_R1_001.fastq.gz    CHRF_R1_gDNA_R1.fastq.gz
CHRF_R1_gDNA_S4_L001_R2_001.fastq.gz    CHRF_R1_gDNA_S4_L002_R2_001.fastq.gz    CHRF_R1_gDNA_R2.fastq.gz
CHRF_R1_gDNA_S4_L001_R3_001.fastq.gz    CHRF_R1_gDNA_S4_L002_R3_001.fastq.gz    CHRF_R1_gDNA_R3.fastq.gz
CHRF_R1_gDNA_S4_L001_I1_001.fastq.gz    CHRF_R1_gDNA_S4_L002_I1_001.fastq.gz    CHRF_R1_gDNA_I1.fastq.gz
CHRF_R1_cDNA_S9_L001_R1_001.fastq.gz    CHRF_R1_cDNA_S9_L002_R1_001.fastq.gz    CHRF_R1_cDNA_R1.fastq.gz
CHRF_R1_cDNA_S9_L001_R2_001.fastq.gz    CHRF_R1_cDNA_S9_L002_R2_001.fastq.gz    CHRF_R1_cDNA_R2.fastq.gz

- The folder with all the original fastq files is in /nfs/team151_data03/MPRA_Programmed/

## SOFTWARE DEPENDENCIES FOR SAMPLE PROCESSING

-Trimmomatic-0.39
-umi_tools version: 1.1.1
-flash2 v2.2.00
-bwa
-samtools
-bamtools 2.5.1


# The barcode reference to align the fastq "Library_TRIMMED_15bp_Carried_Variants.fasta" is in the Dependencies folder with the indexes needed by bwa.
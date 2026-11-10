#!/usr/bin/env Rscript
## Usage:
# Rscript perform_bcalm_variants_cli.R \
#     ../cCRE_based_analysis/data/reporter_experiment.barcode.NGN2.bbmapMapq30StrandSensitiveAssignment.default.all.tsv.gz \
#     data/2509_variant_map_region_based_controls.tsv.gz \
#     <Path to BCalm> \
#     3 \
#     data/ \
#     NGN2_variant_bcalm_80k_bbmap_251128_bcalm_normalized_variant_map_names_dev.tsv \
#     TRUE \
#     10 \
#     1 \
#     2

# NOTE: You have to clone the BCalm repository locally and put the absolute path to the BCalm directory in the commmand above (https://github.com/kircherlab/BCalm) and if you want to use the same version 0.9.0 was used.
# git clone git@github.com:kircherlab/BCalm.git

## Note if you need to install packages, uncomment:
# install.packages("dplyr")
# install.packages("devtools")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager", repos = 'https://cran.uni-muenster.de/')
# BiocManager::install("mpra")
# install.packages("dplyr")
# install.packages("devtools")
# install.packages("tidyr")
# # install ggplot
# install.packages("ggplot2")

# -----------------------------
# Load packages
# -----------------------------
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(devtools)
    library(tidyr)
    library(arrow)
})

# -----------------------------
# Parse command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 10) {
    stop("Usage: Rscript perform_bcalm_variants_kisa_server.R <input> <variant_map> <bcalm_path> <n_reps> <output_dir> <output_name> <normalize> <bc_per_oligo_threshold> <bc_dna_threshold> <bc_found_in_n_reps>")
}

input <- args[1]
variant_map <- args[2]
bcalm_path <- args[3]
n_reps <- as.integer(args[4])
output_dir <- args[5]
output_name <- args[6]
normalize <- ifelse(length(args) >= 7, as.logical(args[7]), TRUE)
bc_per_oligo_threshold <- ifelse(length(args) >= 8, as.integer(args[8]), 10)
bc_dna_threshold <- ifelse(length(args) >= 9, as.integer(args[9]), 1)
bc_found_in_n_reps <- ifelse(length(args) >= 10, as.integer(args[10]), 2)

# -----------------------------
# Load BCalm package
# -----------------------------
load_all(bcalm_path)

# -----------------------------
# Read input
# -----------------------------
message("Reading input file: ", input)
df <- read.table(file = gzfile(input), sep = "\t", header = TRUE, comment.char = "")

colnames(df) <- c(
    "barcode", "name", "dna_count_1", "rna_count_1",
    "dna_count_2", "rna_count_2", "dna_count_3", "rna_count_3"
)


# -----------------------------
# Filter for positive counts
# -----------------------------
df <- df %>%
    filter(
        rowSums(across(matches("DNA", ignore.case = TRUE), ~ . >= bc_dna_threshold)) >= bc_found_in_n_reps,
    )

message("Number of unique names after filtering: ", length(unique(df$name)))

# -----------------------------
# Variant map
# -----------------------------
message("Reading variant map: ", variant_map)
map_df <- read.table(file = gzfile(variant_map), sep = "\t", header = TRUE, comment.char = "")
map_df <- map_df %>% select(ID, REF, ALT)

# -----------------------------
# Merge for REF and ALT
# -----------------------------
df_ref <- merge(df, map_df, by.x = "name", by.y = "REF", all.x = FALSE)
df_ref$allele <- "ref"
df_ref$ALT <- NULL
message("REF oligos: ", length(unique(df_ref$name)))

df_alt <- merge(df, map_df, by.x = "name", by.y = "ALT", all.x = FALSE)
df_alt$allele <- "alt"
df_alt$REF <- NULL
message("ALT oligos: ", length(unique(df_alt$name)))

df_combined <- rbind(df_ref, df_alt)

var_df <- df_combined %>% select(variant_id = ID, allele, barcode, matches("count"))
message("Unique variant IDs: ", length(unique(var_df$variant_id)))

# -----------------------------
# Filter for min barcodes and both alleles
# -----------------------------
df_filt <- var_df %>%
    group_by(variant_id, allele) %>%
    filter(n() >= bc_per_oligo_threshold) %>%
    ungroup() %>%
    group_by(variant_id) %>%
    filter(n_distinct(allele) == 2) %>%
    ungroup()


message("Variants after barcode/allele filtering: ", length(unique(df_filt$variant_id)))

# -----------------------------
# Downsample and create MPRA input
# -----------------------------
df_smaller <- downsample_barcodes(df_filt, id_column_name = "variant_id")
dna <- create_dna_df(df_smaller, id_column_name = "variant_id")
rna <- create_rna_df(df_smaller, id_column_name = "variant_id")

bcs <- ncol(dna) / (n_reps * 2)
# print Barcodes per oligo and n_reps
message("Barcodes per oligo: ", bcs, ", Number of replicates: ", n_reps)

# -----------------------------
# Run BCalm
# -----------------------------
message("Creating MPRASet and running mpralm...")

mpra <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL)
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpra)))
block_vector <- rep(1:n_reps, each = bcs * 2)

start <- Sys.time()
mpralm_fit <- mpralm(
    object = mpra, design = design, aggregate = "none",
    normalize = normalize, model_type = "corr_groups", plot = FALSE,
    block = block_vector
)
message("Finished in ", round(Sys.time() - start, 2), " seconds")

# -----------------------------
# Write output
# -----------------------------
output_path <- file.path(output_dir, output_name)
message("Writing results to: ", output_path)

toptab_allele_bc <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab_allele_bc$variant_id <- row.names(toptab_allele_bc)
write.table(toptab_allele_bc, file = output_path, sep = "\t", row.names = FALSE)

message("Done")

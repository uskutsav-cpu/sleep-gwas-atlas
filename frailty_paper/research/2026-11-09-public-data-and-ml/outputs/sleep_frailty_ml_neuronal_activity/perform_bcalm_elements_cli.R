#!/usr/bin/env Rscript

## Usage:
# Rscript perform_bcalm_elements_cli.R \
#     <Path to BCalm> \
#     data/reporter_experiment.barcode.NGN2.bbmapMapq30StrandSensitiveAssignment.default.all.tsv.gz \
#     data/mpra_80k_seq_label.tsv.gz \
#     0.95 \
#     negative_neuron_ctrl \
#     scramble_ctrl \
#     cardiac_neuro_cava_random \
#     data/cli_NGN2_202511_80k_element_bbmap_bcalm_normalized_names_bc10_toptreat_results_negative_ctrl_result_095_dev.tsv \
#     data/cli_NGN2_202511_80k_element_bbmap_bcalm_normalized_bc10_toptreat_results_scrambled_NP_MK_result_095_dev.tsv \
#     data/cli_NGN2_202511_80k_element_bbmap_bcalm_normalized_names_negative_ctrl_ctrl_095.png \
#     data/cli_NGN2_202511_80k_element_bbmap_bcalm_normalized_names_scrambled_ctrl_095.png \
#     3 \
#     TRUE \
#     10 \
#     1 \
#     2


# NOTE: You have to clone the BCalm repository locally and put the absolute path to the BCalm directory in the commmand above (https://github.com/kircherlab/BCalm) and if you want to use the same version 0.9.0 was used.
# git clone git@github.com:kircherlab/BCalm.git


## NOTE: if you need to install packages uncomment:
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


# -------------------------------------
# Load libraries
# -------------------------------------
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(devtools)
    library(tidyr)
})

# -------------------------------------
# Helper: timestamped log messages
# -------------------------------------
log_msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "-", ..., "\n")

# -------------------------------------
# Parse command-line arguments
# -------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 12) {
    cat("
Usage:
  Rscript run_bcalm_element_analysis.R \\
    <bcalm_path> <input_counts> <label_file> <sig_level> \\
    <negative_ctrl_label> <scramble_ctrl_label> <test_label> \\
    <toptreat_negative_out> <toptreat_scrambled_out> \\
    <plot_negative_out> <plot_scrambled_out> <nr_reps>

Arguments:
  <bcalm_path>             Path to the local BCalm installation directory
  <input_counts>           Path to normalized element count table (.tsv.gz)
  <label_file>             Path to label file (must contain 'name' and 'label' columns)
  <sig_level>              Significance level (e.g., 0.95)
  <negative_ctrl_label>    Label of the negative control sequences
  <scramble_ctrl_label>    Label of the scrambled control sequences (if available)
  <test_label>             Label of the test group to compare
  <toptreat_negative_out>  Output path for topTreat results vs negative control (.tsv)
  <toptreat_scrambled_out> Output path for topTreat results vs scrambled control (.tsv)
  <plot_negative_out>      Output path for plot vs negative control (.png)
  <plot_scrambled_out>     Output path for plot vs scrambled control (.png)
  <nr_reps>                Number of replicates (integer)
  <normalize>              Whether to normalize data (TRUE/FALSE)
  <bc_per_oligo_threshold> Oligo DNA threshold for filtering low count oligos out (default: 10)
  <bc_dna_threshold>       DNA threshold for filtering barcodes with low counts (default: 1)
  <bc_found_in_n_reps>     Minimum number of replicates a barcode must be found in (default: 2)

Usage:
  Rscript run_bcalm_element_analysis.R \\
    /home/user/projects/BCalm \\
    data/80k_element_counts.tsv.gz \\
    data/labels.tsv.gz \\
    0.95 \\
    negative_ctrl scramble_ctrl cardiac_neuro_cava_random \\
    results/toptreat_vs_neg.tsv results/toptreat_vs_scrambled.tsv \\
    results/plot_vs_neg.png results/plot_vs_scrambled.png \\
    3 \\
    TRUE \\
    10 \\
    1
\n")
    quit(status = 1)
}

# -------------------------------------
# Assign arguments
# -------------------------------------
bcalm_path <- args[1]
input_counts <- args[2]
labelfile <- args[3]
sig_level_mpra_treat <- as.numeric(args[4])
negative_control <- args[5]
scramble_control <- args[6]
test_label <- args[7]
toptreat_negative_out <- args[8]
toptreat_scrambled_out <- args[9]
plot_negative_out <- args[10]
plot_scrambled_out <- args[11]
nr_reps <- as.integer(args[12])
normalize <- tolower(args[13]) %in% c("true", "t", "1", "yes", "y")
bc_per_oligo_threshold <- ifelse(length(args) >= 14, as.integer(args[14]), 10)
bc_dna_threshold <- ifelse(length(args) >= 15, as.integer(args[15]), 1)
bc_found_in_n_reps <- ifelse(length(args) >= 16, as.integer(args[16]), 2)

# -------------------------------------
# Check input files
# -------------------------------------
if (!dir.exists(bcalm_path)) stop("BCalm path not found: ", bcalm_path)
if (!file.exists(input_counts)) stop("Input count file not found: ", input_counts)
if (!file.exists(labelfile)) stop("Label file not found: ", labelfile)

# -------------------------------------
# Load BCalm package dynamically
# -------------------------------------
log_msg("Loading BCalm package from:", bcalm_path)
load_all(bcalm_path)

# -------------------------------------
# Read normalized element count table
# -------------------------------------
log_msg("Reading label file:", labelfile)

log_msg("Reading count table...", input_counts)
df <- read.table(file = gzfile(input_counts), sep = "\t", header = TRUE, comment.char = "")

colnames(df) <- c(
    "barcode", "name", "dna_count_1", "rna_count_1",
    "dna_count_2", "rna_count_2", "dna_count_3", "rna_count_3"
)

log_msg("Loaded", nrow(df), "rows with", length(unique(df$name)), "unique elements.")

# Filter: oligos with at least 2/3 DNA and RNA replicates > 0
df <- df %>%
    filter(
        rowSums(across(matches("DNA", ignore.case = TRUE), ~ . >= bc_dna_threshold)) >= bc_found_in_n_reps,
        # Do not filter on RNA: changes the overall rna size
        # rowSums(across(matches("RNA", ignore.case = TRUE), ~ . > 0)) >= bc_found_in_n_reps
    )
log_msg("After filtering:", nrow(df), "rows,", length(unique(df$name)), "unique names remain.")

# Barcode threshold filter
df <- df %>%
    group_by(name) %>%
    filter(n() >= bc_per_oligo_threshold) %>%
    ungroup()
log_msg("After barcode filter:", length(unique(df$name)), "unique elements.")

# -------------------------------------
# Create MPRA input data
# -------------------------------------
df_small <- downsample_barcodes(df)
dna <- create_dna_df(df_small, id_column_name = "name")
rna <- create_rna_df(df_small, id_column_name = "name")

bcs <- ncol(dna) / nr_reps
block_vector <- rep(1:nr_reps, each = bcs)

# -------------------------------------
# Read label file (check for header, require name + label)
# -------------------------------------
log_msg("Reading label file:", labelfile)

label_preview <- readLines(gzfile(labelfile), n = 1)
has_header <- grepl("name", label_preview, ignore.case = TRUE) &&
    grepl("label", label_preview, ignore.case = TRUE)

if (has_header) {
    labels <- read.table(gzfile(labelfile), header = TRUE, sep = "\t", comment.char = "")
} else {
    log_msg("No header detected in label file; assuming first two columns are 'name' and 'label'.")
    labels <- read.table(gzfile(labelfile), header = FALSE, sep = "\t", comment.char = "")
    if (ncol(labels) < 2) stop("Label file must have at least two columns: 'name' and 'label'.")
    colnames(labels)[1:2] <- c("name", "label")
}

required_cols <- c("name", "label")
if (!all(required_cols %in% colnames(labels))) {
    stop("Label file must contain columns: ", paste(required_cols, collapse = ", "))
}

log_msg("Label file columns:", paste(colnames(labels), collapse = ", "))

labels_vec <- setNames(labels$label, labels$name)
labels_vec <- labels_vec[rownames(dna)]
log_msg("Matched", sum(!is.na(labels_vec)), "labels to DNA matrix rows.")

# -------------------------------------
# Run BCalm model
# -------------------------------------
log_msg("Running fit_elements (BCalm model)...")
mpra <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL, label = labels_vec)
start <- Sys.time()
mpralm_fit <- fit_elements(object = mpra, normalize = normalize, block = block_vector, plot = FALSE)
log_msg("BCalm fit complete in", round(Sys.time() - start, 2), "seconds.")

# -------------------------------------
# Compare vs negative control
# -------------------------------------
log_msg("Running mpra_treat for negative control:", negative_control)
testing_element_negative_control <- mpra_treat(mpralm_fit, sig_level_mpra_treat, neg_label = negative_control)
element_toptreat_neg <- topTreat(testing_element_negative_control, coef = 1, number = Inf)
element_toptreat_neg$name <- row.names(element_toptreat_neg)

log_msg("Writing topTreat results (negative control):", toptreat_negative_out)
write.table(element_toptreat_neg,
    file = toptreat_negative_out,
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
)

log_msg("Generating plot for negative control comparison...")
element_plot_neg <- plot_groups(mpralm_fit, sig_level_mpra_treat, neg_label = negative_control, test_label = test_label)
ggsave(filename = plot_negative_out, plot = element_plot_neg, width = 6, height = 5, dpi = 300)
log_msg("Saved plot:", plot_negative_out)

# -------------------------------------
# Compare vs scrambled control
# -------------------------------------
log_msg("Running mpra_treat for scrambled control:", scramble_control)
testing_element_scrambled <- mpra_treat(mpralm_fit, sig_level_mpra_treat, neg_label = scramble_control)
element_toptreat_scr <- topTreat(testing_element_scrambled, coef = 1, number = Inf)
element_toptreat_scr$name <- row.names(element_toptreat_scr)

log_msg("Writing topTreat results (scrambled control):", toptreat_scrambled_out)
write.table(element_toptreat_scr,
    file = toptreat_scrambled_out,
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
)

log_msg("Generating plot for scrambled control comparison...")
element_plot_scr <- plot_groups(mpralm_fit, sig_level_mpra_treat, neg_label = scramble_control, test_label = test_label)
ggsave(filename = plot_scrambled_out, plot = element_plot_scr, width = 6, height = 5, dpi = 300)
log_msg("Saved plot:", plot_scrambled_out)

# -------------------------------------
# Finish
# -------------------------------------
log_msg("Analysis complete.")
log_msg("Results written to:")
log_msg("-", toptreat_negative_out)
log_msg("-", toptreat_scrambled_out)
log_msg("Plots saved to:")
log_msg("-", plot_negative_out)
log_msg("-", plot_scrambled_out)

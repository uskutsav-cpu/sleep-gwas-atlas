# Load libraries
library(tidyverse)
library(data.table)
library(BuenColors)
library(binom)
library(ggplot2)
library(cowplot)
library(mashr)
library(patchwork)
library("PNWColors")
source("code/utils.R")
options(stringsAsFactors = FALSE)

# read in processed MPRA data
mpra_df <- vroom('data/preprocess/core_mpra.txt.gz') %>% 
  dplyr::filter(cohort == 'control' | (pip > 0.1 | cs_id > 0))

# Get moderate-high PIP variants
df <- mpra_df %>%
  dplyr::filter(cohort %in% c("UKBB", "GTEx", "BBJ")) %>%
  dplyr::mutate(cohort = ifelse(cohort %in% c("UKBB", "BBJ"), "complex", cohort)) %>%
  dplyr::filter(CRE > 0, pip >= 0.1) %>%
  dplyr::group_by(cohort, variant, cell_type) %>%
  dplyr::filter(pip == max(pip),
                !is.na(pip)) %>%
  dplyr::filter(row_number() == 1) %>%
  ungroup() 

# Create beta / se matrices
ukbb_beta_df <- df %>%
  dplyr::filter(cohort %in% "complex") %>%
  dplyr::select(variant, cell_type, log2Skew) %>%
  pivot_wider(names_from = "cell_type",
              values_from = log2Skew) %>%
  drop_na()
gtex_beta_df <- df %>%
  dplyr::filter(cohort %in% c("GTEx")) %>%
  dplyr::select(variant, cell_type, log2Skew) %>%
  pivot_wider(names_from = "cell_type",
              values_from = log2Skew) %>%
  drop_na()

ukbb_se_df <- df %>%
  dplyr::filter(cohort %in% c("complex")) %>%
  dplyr::select(variant, cell_type, Skew_SE) %>%
  pivot_wider(names_from = "cell_type",
              values_from = Skew_SE) %>%
  drop_na()
gtex_se_df <- df %>%
  dplyr::filter(cohort %in% c("GTEx")) %>%
  dplyr::select(variant, cell_type, Skew_SE) %>%
  pivot_wider(names_from = "cell_type",
              values_from = Skew_SE) %>%
  drop_na()

ukbb_skew_df <- df %>%
  dplyr::filter(cohort %in% c("complex")) %>%
  dplyr::select(variant, cell_type, emVar) %>%
  pivot_wider(names_from = "cell_type",
              values_from = emVar) %>%
  drop_na() %>%
  dplyr::mutate(skew = ifelse(K562 | HEPG2 | SKNSH | A549, T, F))
gtex_skew_df <- df %>%
  dplyr::filter(cohort %in% c("GTEx")) %>%
  dplyr::select(variant, cell_type, emVar) %>%
  pivot_wider(names_from = "cell_type",
              values_from = emVar) %>%
  dplyr::mutate(skew = ifelse(HEPG2 | SKNSH | HCT116, T, F)) %>%
  drop_na()

ukbb_data <- mash_set_data(as.matrix(ukbb_beta_df[,-1] %>% data.table(), rownames.value = ukbb_beta_df$variant),
                           as.matrix(ukbb_se_df[,-1] %>% data.table(), rownames.value = ukbb_se_df$variant))
gtex_data <- mash_set_data(as.matrix(gtex_beta_df[,-1] %>% data.table(), rownames.value = gtex_beta_df$variant),
                           as.matrix(gtex_se_df[,-1] %>% data.table(), rownames.value = gtex_se_df$variant))
# Fit mash model
ukbb_U.c <- cov_canonical(ukbb_data)  
gtex_U.c <- cov_canonical(gtex_data)  
ukbb_m.c <- mash(ukbb_data, ukbb_U.c)
gtex_m.c <- mash(gtex_data, gtex_U.c)

# Get model outputs
ukbb_pm <- ukbb_m.c$result$PosteriorMean %>%
  data.frame() %>%
  setNames(paste0('pm.', names(.))) %>%
  rownames_to_column("variant") %>%
  as_tibble() 
gtex_pm <- gtex_m.c$result$PosteriorMean %>%
  data.frame() %>%
  setNames(paste0('pm.', names(.))) %>%
  rownames_to_column("variant") %>%
  as_tibble() 

ukbb_pd <- ukbb_m.c$result$PosteriorSD %>%
  data.frame() %>%
  setNames(paste0('sd.', names(.))) %>%
  as_tibble() 
gtex_pd <- gtex_m.c$result$PosteriorSD %>%
  data.frame() %>%
  setNames(paste0('sd.', names(.))) %>%
  as_tibble() 

ukbb_lfsr <- ukbb_m.c$result$lfsr %>%
  data.frame() %>%
  setNames(paste0('lfsr.', names(.))) %>%
  as_tibble() 
gtex_lfsr <- gtex_m.c$result$lfsr %>%
  data.frame() %>%
  setNames(paste0('lfsr.', names(.))) %>%
  as_tibble() 

# Combine outputs into tidy data frame
ukbb_mash_df <- bind_cols(ukbb_pm, ukbb_pd, ukbb_lfsr) %>%
  pivot_longer(cols = !variant) %>%
  separate(name, into = c("estimate", "cell_type"), sep = "\\.") %>%
  pivot_wider(names_from = "estimate", values_from = "value")
gtex_mash_df <- bind_cols(gtex_pm, gtex_pd, gtex_lfsr) %>%
  pivot_longer(cols = !variant) %>%
  separate(name, into = c("estimate", "cell_type"), sep = "\\.") %>%
  pivot_wider(names_from = "estimate", values_from = "value")

# Get cell-type specific and non-specific emVars effects
# lfsr > 1 is non-specific, lfsr == 1 is specific
ukbb_mash_df <- ukbb_mash_df %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(num_emVar = sum(lfsr < 0.05))
gtex_mash_df <- gtex_mash_df %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(num_emVar = sum(lfsr < 0.05))

rbind(ukbb_mash_df %>% dplyr::mutate(cohort = 'Complex Traits'), 
      gtex_mash_df %>% dplyr::mutate(cohort = 'eQTLs')) %>%
  vroom::vroom_write("tables/stable7.skew.txt.gz", delim = '\t', col_names = T)


table(ukbb_mash_df$num_emVar)[-1] / sum(table(ukbb_mash_df$num_emVar)[-1])
table(gtex_mash_df$num_emVar)[-1] / sum(table(gtex_mash_df$num_emVar)[-1])


# script to run parameter sensitivity analysis on deMULTIplex2

##########################################
#      TUNABLE PARAMS:
#        *init.cos.cut — initial cosine similarity cutoff (default 0.5)
#        *prob.cut — probability cutoff for classification (default 0.5)
#        min.cell.fit — minimum cells to fit model (default 10)
#        converge.threshold — convergence threshold for EM (default 0.001)
##########################################

# TO RUN:
#        conda activate demux-r
#        Rscript analysis/parameter_analysis/sensitivity_dm2.R

library(deMULTIplex2)
library(tidyverse)

# load scoring functions
source('analysis/parameter_analysis/scores.R')

# load ground truth and scoring inputs 
load('analysis/parameter_analysis/mcginnis_ms_scores.RData')
load('analysis/parameter_analysis/mcginnis_ab_scores.RData')
truth_ms <- truth_mcginnis_ms_vireo
bars_ms <- bars_mcginnis_ms
weights_ms <- weights_mcginnis_ms

truth_ab <- truth_mcginnis_ab_vireo
bars_ab <- bars_mcginnis_ab
weights_ab <- weights_mcginnis_ab

# load HTO data
data_ms <- read.csv('data/mcginnis_ms/hto/GSM4904942_8donor_PBMC_AH_MULTI_matrix.csv', row.names = 1)
data_ms <- data_ms[, !colnames(data_ms) %in% c('nUMI', 'nUMI_total')]
data_ms <- as.matrix(data_ms[, 1:8])

data_ab <- read.csv('data/mcginnis_hto/hto/GSM4904939_8donor_PBMC_AH_SCMK_matrix.csv', row.names = 1)
data_ab <- data_ab[, !colnames(data_ab) %in% c('nUMI', 'nUMI_total')]
data_ab <- as.matrix(data_ab[, 1:8])

# load barcode maps
barcode_map_ms <- read.csv('config/barcode_maps/mcginnis_ms.csv')
barcode_lookup_ms <- setNames(barcode_map_ms$index, barcode_map_ms$barcode)

barcode_map_ab <- read.csv('config/barcode_maps/mcginnis_ab.csv')
barcode_lookup_ab <- setNames(barcode_map_ab$index, barcode_map_ab$barcode)

# parameter grid
init_cos_vals <- c(0.3, 0.4, 0.5, 0.6, 0.7)
prob_cut_vals <- c(0.3, 0.4, 0.5, 0.6, 0.7)

#parameter scan MCGINNIS_MS
dm2_results_ms <- list()
i <- 1

for (init_cos in init_cos_vals) {
  for (prob_cut in prob_cut_vals) {
    
    res <- tryCatch({
      demultiplexTags(data_ms,
                      init.cos.cut = init_cos,
                      prob.cut = prob_cut,
                      plot.umap = 'none',
                      plot.diagnostics = FALSE)
    }, error = function(e) {
      message(paste('Failed ms:', init_cos, prob_cut, e$message))
      NULL
    })
    
    if (is.null(res)) next
    
    pred_named <- res$final_assign
    pred <- case_when(
      str_detect(pred_named, 'negative') ~ 0,
      str_detect(pred_named, 'multiplet') ~ 1000,
      pred_named %in% names(barcode_lookup_ms) ~ barcode_lookup_ms[pred_named],
      TRUE ~ NA_real_
    )
    names(pred) <- rownames(res$assign_table)
    
    common_cells <- intersect(names(truth_ms), names(pred))
    tr <- truth_ms[common_cells]
    pr <- pred[common_cells]
    
    dm2_results_ms[[i]] <- data.frame(
      init_cos_cut = init_cos,
      prob_cut = prob_cut,
      f1_micro = f1_micro(tr, pr),
      f1_macro = f1_macro(tr, pr, bars_ms),
      f1_weighted = f1_weighted(tr, pr, bars_ms, weights_ms),
      mcc = mcc(tr, pr, bars_ms),
      concordance = concordance(tr, pr)
    )
    i <- i + 1
    message(paste('Done ms:', init_cos, prob_cut))
  }
}

dm2_results_ms_df <- bind_rows(dm2_results_ms)
write.csv(dm2_results_ms_df, 'analysis/parameter_analysis/dm2_mcginnis_ms.csv', row.names = FALSE)

#parameter scan MCGINNIS_AB
dm2_results_ab <- list()
i <- 1

for (init_cos in init_cos_vals) {
  for (prob_cut in prob_cut_vals) {
    
    res <- tryCatch({
      demultiplexTags(data_ab,
                      init.cos.cut = init_cos,
                      prob.cut = prob_cut,
                      plot.umap = 'none',
                      plot.diagnostics = FALSE)
    }, error = function(e) {
      message(paste('Failed ab:', init_cos, prob_cut, e$message))
      NULL
    })

    if (is.null(res)) next
    
    pred_named <- res$final_assign
    pred <- case_when(
      str_detect(pred_named, 'negative') ~ 0,
      str_detect(pred_named, 'multiplet') ~ 1000,
      pred_named %in% names(barcode_lookup_ab) ~ barcode_lookup_ab[pred_named],
      TRUE ~ NA_real_
    )
    names(pred) <- rownames(res$assign_table)

    common_cells <- intersect(names(truth_ab), names(pred))
    tr <- truth_ab[common_cells]
    pr <- pred[common_cells]

    dm2_results_ab[[i]] <- data.frame(
      init_cos_cut = init_cos,
      prob_cut = prob_cut,
      f1_micro = f1_micro(tr, pr),
      f1_macro = f1_macro(tr, pr, bars_ab),
      f1_weighted = f1_weighted(tr, pr, bars_ab, weights_ab),
      mcc = mcc(tr, pr, bars_ab),
      concordance = concordance(tr, pr)
    )
    i <- i + 1
    message(paste('Done ab:', init_cos, prob_cut))
  }
}

dm2_results_ab_df <- bind_rows(dm2_results_ab)
write.csv(dm2_results_ab_df, 'analysis/parameter_analysis/dm2_mcginnis_ab.csv', row.names = FALSE)
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
source('analysis/parameter_analysis/scores.R')  # or paste the functions directly

# load ground truth and scoring inputs
load('analysis/parameter_analysis/mcginnis_ms_scores.RData')
truth <- truth_mcginnis_ms_vireo
bars <- bars_mcginnis_ms
weights <- weights_mcginnis_ms

# load HTO data
data <- read.csv('data/mcginnis_ms/hto/GSM4904942_8donor_PBMC_AH_MULTI_matrix.csv', row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')]
data <- as.matrix(data[, 1:8])

# load barcode map
barcode_map <- read.csv('config/barcode_maps/mcginnis_ms.csv')
barcode_lookup <- setNames(barcode_map$index, barcode_map$barcode)

#parameter scan
init_cos_vals <- c(0.3, 0.4, 0.5, 0.6, 0.7)
prob_cut_vals <- c(0.3, 0.4, 0.5, 0.6, 0.7)

dm2_results <- list()
i <- 1

for (init_cos in init_cos_vals) {
  for (prob_cut in prob_cut_vals) {
    
    res <- tryCatch({
      demultiplexTags(data,
                      init.cos.cut = init_cos,
                      prob.cut = prob_cut,
                      plot.umap = 'none',
                      plot.diagnostics = FALSE)
    }, error = function(e) {
      message(paste('Failed:', init_cos, prob_cut, e$message))
      NULL
    })
    
    if (is.null(res)) next
    
    pred_named <- res$final_assign
    pred <- case_when(
      str_detect(pred_named, 'negative') ~ 0,
      str_detect(pred_named, 'multiplet') ~ 1000,
      pred_named %in% names(barcode_lookup) ~ barcode_lookup[pred_named],
      TRUE ~ NA_real_
    )
    names(pred) <- rownames(res$assign_table)
    
    common_cells <- intersect(names(truth), names(pred))
    tr <- truth[common_cells]
    pr <- pred[common_cells]
    
    dm2_results[[i]] <- data.frame(
      init_cos_cut = init_cos,
      prob_cut = prob_cut,
      f1_micro = f1_micro(tr, pr),
      f1_macro = f1_macro(tr, pr, bars),
      f1_weighted = f1_weighted(tr, pr, bars, weights),
      mcc = mcc(tr, pr, bars),
      concordance = concordance(tr, pr)
    )
    i <- i + 1
    message(paste('Done:', init_cos, prob_cut))
  }
}

dm2_results_df <- bind_rows(dm2_results)
write.csv(dm2_results_df, 'analysis/parameter_analysis/dm2_mcginnis_ms.csv', row.names = FALSE)
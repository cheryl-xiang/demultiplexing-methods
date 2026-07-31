#Script to calculate scores for HashSolo parameter scan

# TO RUN:
#        conda activate demux-r
#        Rscript analysis/parameter_analysis/hsolo_scoring.R

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

score_runs <- function(files, truth, bars, weights) {
  results <- list()
  i <- 1
  for (f in files) {
    run <- read.csv(f)
    pred <- setNames(run$assignment, run$cell_barcode)
    
    common_cells <- intersect(names(truth), names(pred))
    tr <- truth[common_cells]
    pr <- pred[common_cells]
    
    results[[i]] <- data.frame(
      prior_neg = run$prior_neg[1],
      prior_singlet = run$prior_singlet[1],
      prior_doublet = run$prior_doublet[1],
      n_noise = run$n_noise[1],
      f1_micro = f1_micro(tr, pr),
      f1_macro = f1_macro(tr, pr, bars),
      f1_weighted = f1_weighted(tr, pr, bars, weights),
      mcc = mcc(tr, pr, bars),
      concordance = concordance(tr, pr)
    )
    i <- i + 1
  }
  bind_rows(results)
}

# ---- mcginnis_ms ----
files_ms <- list.files('analysis/parameter_analysis/hashsolo_runs/mcginnis_ms', full.names = TRUE)
hashsolo_ms_df <- score_runs(files_ms, truth_ms, bars_ms, weights_ms)
write.csv(hashsolo_ms_df, 'analysis/parameter_analysis/hashsolo_mcginnis_ms.csv', row.names = FALSE)
message('Done mcginnis_ms')

# ---- mcginnis_ab ----
files_ab <- list.files('analysis/parameter_analysis/hashsolo_runs/mcginnis_ab', full.names = TRUE)
hashsolo_ab_df <- score_runs(files_ab, truth_ab, bars_ab, weights_ab)
write.csv(hashsolo_ab_df, 'analysis/parameter_analysis/hashsolo_mcginnis_ab.csv', row.names = FALSE)
message('Done mcginnis_ab')

message('Done :)')
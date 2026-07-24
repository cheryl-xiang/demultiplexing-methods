#Script to calculate scores for HashSolo parameter scan

# TO RUN:
#        conda activate demux-r
#        Rscript analysis/parameter_analysis/hsolo_scoring.R

library(tidyverse)

# load scoring functions
source('analysis/parameter_analysis/scores.R')  # or paste the functions directly

# load ground truth and scoring inputs
load('analysis/parameter_analysis/mcginnis_ms_scores.RData')
truth <- truth_mcginnis_ms_vireo
bars <- bars_mcginnis_ms
weights <- weights_mcginnis_ms

# load all hashsolo runs
files <- list.files('analysis/parameter_analysis/hashsolo_runs', full.names = TRUE)

hashsolo_results <- list()
i <- 1

for (f in files) {
  run <- read.csv(f)
  pred <- setNames(run$assignment, run$cell_barcode)
  
  common_cells <- intersect(names(truth), names(pred))
  tr <- truth[common_cells]
  pr <- pred[common_cells]
  
  hashsolo_results[[i]] <- data.frame(
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

hashsolo_df <- bind_rows(hashsolo_results)
write.csv(hashsolo_df, 'analysis/parameter_analysis/hashsolo_mcginnis_ms.csv', row.names = FALSE)

message(paste('Done :)'))
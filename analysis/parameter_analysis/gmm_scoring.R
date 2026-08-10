#Script to calculate scores for GMMDemux calls against ground truth

# TO RUN:
#        conda activate demux-r
#        Rscript analysis/parameter_analysis/gmm_scoring.R

library(tidyverse)

source('analysis/parameter_analysis/scores.R')

load('analysis/parameter_analysis/scores/mcginnis_ms_scores.RData')
load('analysis/parameter_analysis/scores/mcginnis_ab_scores.RData')

truth_ms <- truth_mcginnis_ms_vireo
bars_ms <- bars_mcginnis_ms
weights_ms <- weights_mcginnis_ms

truth_ab <- truth_mcginnis_ab_vireo
bars_ab <- bars_mcginnis_ab
weights_ab <- weights_mcginnis_ab

score_assignments <- function(assignments_file, truth, bars, weights, dataset, method) {
  if (!file.exists(assignments_file)) {
    message(paste('File not found:', assignments_file))
    return(NULL)
  }
  
  assignments <- read.csv(assignments_file)
  pred <- setNames(assignments$assignment, assignments$cell_barcode)
  
  common_cells <- intersect(names(truth), names(pred))
  tr <- truth[common_cells]
  pr <- pred[common_cells]
  
  data.frame(
    dataset = dataset,
    method = method,
    f1_micro = f1_micro(tr, pr),
    f1_macro = f1_macro(tr, pr, bars),
    f1_weighted = f1_weighted(tr, pr, bars, weights),
    mcc = mcc(tr, pr, bars),
    concordance = concordance(tr, pr)
  )
}

results <- bind_rows(
  score_assignments('results/gmmdemux/mcginnis_ms/assignments.csv',
                    truth_ms, bars_ms, weights_ms, 'mcginnis_ms', 'gmmdemux'),
  score_assignments('results/gmmdemux/mcginnis_hto/assignments.csv',
                    truth_ab, bars_ab, weights_ab, 'mcginnis_ab', 'gmmdemux')
)

print(results)
write.csv(results, 'analysis/parameter_analysis/gmmdemux_scores.csv', row.names = FALSE)
message('Done :)')
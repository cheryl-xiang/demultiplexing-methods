#Script to calculate scores for Beacon calls against ground truth

# TO RUN:
#        conda activate demux-r
#        Rscript analysis/beacon_analysis/beacon_scoring.R

library(tidyverse)

source('analysis/parameter_analysis/scores.R')

#load scores
load('analysis/parameter_analysis/scores/mcginnis_ms_scores.RData')
load('analysis/parameter_analysis/scores/mcginnis_ab_scores.RData')
load('analysis/parameter_analysis/scores/bar11_scores.RData')
load('analysis/parameter_analysis/scores/gaublomme_scores.RData')
load('analysis/parameter_analysis/scores/howitt_cell_line_scores.RData')
load('analysis/parameter_analysis/scores/howitt_scores.RData')

#get truths, bars, weights
truth_ms <- truth_mcginnis_ms_vireo
bars_ms <- bars_mcginnis_ms
weights_ms <- weights_mcginnis_ms

truth_ab <- truth_mcginnis_ab_vireo
bars_ab <- bars_mcginnis_ab
weights_ab <- weights_mcginnis_ab

truth_bar11 <- truth_bar11_vireo

truth_gaublomme <- truth_gaublomme_demuxlet

truth_howitt_cap1 <- truth_howitt_capture1_vireo
bars_howitt_cap1 <- bars_howitt_capture1
weights_howitt_cap1 <- weights_howitt_capture1

truth_howitt_cap2 <- truth_howitt_capture2_vireo
bars_howitt_cap2 <- bars_howitt_capture2
weights_howitt_cap2 <- weights_howitt_capture2

truth_howitt_cap3 <- truth_howitt_capture3_vireo
bars_howitt_cap3 <- bars_howitt_capture3
weights_howitt_cap3 <- weights_howitt_capture3

truth_howitt_b1c1 <- truth_howitt_batch1_c1_vireo
bars_howitt_b1c1 <- bars_howitt_batch1_c1
weights_howitt_b1c1 <- weights_howitt_batch1_c1

truth_howitt_b1c2 <- truth_howitt_batch1_c2_vireo
bars_howitt_b1c2 <- bars_howitt_batch1_c2
weights_howitt_b1c2 <- weights_howitt_batch1_c2

truth_howitt_b2c1 <- truth_howitt_batch2_c1_vireo
bars_howitt_b2c1 <- bars_howitt_batch2_c1
weights_howitt_b2c1 <- weights_howitt_batch2_c1

truth_howitt_b2c2 <- truth_howitt_batch2_c2_vireo
bars_howitt_b2c2 <- bars_howitt_batch2_c2
weights_howitt_b2c2 <- weights_howitt_batch2_c2

truth_howitt_b3c1 <- truth_howitt_batch3_c1_vireo
bars_howitt_b3c1 <- bars_howitt_batch3_c1
weights_howitt_b3c1 <- weights_howitt_batch3_c1

truth_howitt_b3c2 <- truth_howitt_batch3_c2_vireo
bars_howitt_b3c2 <- bars_howitt_batch3_c2
weights_howitt_b3c2 <- weights_howitt_batch3_c2

#calculate scores
score_assignments <- function(assignments_file, truth, bars, weights, dataset, method, strip_suffix = NULL) {
  if (!file.exists(assignments_file)) {
    message(paste('File not found:', assignments_file))
    return(NULL)
  }
  
  assignments <- read.csv(assignments_file)
  
  # strip leading X\d+_ prefix and trailing .\d+ suffix
  assignments$cell_barcode <- sub('^X\\d+_', '', assignments$cell_barcode)
  assignments$cell_barcode <- sub('\\.\\d+$', '', assignments$cell_barcode)
  
  if (!is.null(strip_suffix)) {
    assignments$cell_barcode <- sub(paste0('\\', strip_suffix, '$'), '', assignments$cell_barcode)
  }
  
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
  score_assignments('results/beacon/mcginnis_ms/assignments.csv',
                    truth_ms, bars_ms, weights_ms, 'mcginnis_ms', 'beacon'),
  score_assignments('results/beacon/mcginnis_hto/assignments.csv',
                    truth_ab, bars_ab, weights_ab, 'mcginnis_ab', 'beacon'),
  score_assignments('results/beacon/bar11/assignments.csv',
                    truth_bar11, bars_bar11, weights_bar11, 'bar11', 'beacon'),
  score_assignments('results/beacon/gaublomme/assignments.csv',
                    truth_gaublomme, bars_gaublomme, weights_gaublomme, 'gaublomme', 'beacon'),
  score_assignments('results/beacon/howitt_capture1/assignments.csv',
                    truth_howitt_cap1, bars_howitt_cap1, weights_howitt_cap1, 'howitt_capture1', 'beacon'),
  score_assignments('results/beacon/howitt_capture2/assignments.csv',
                    truth_howitt_cap2, bars_howitt_cap2, weights_howitt_cap2, 'howitt_capture2', 'beacon'),
  score_assignments('results/beacon/howitt_capture3/assignments.csv',
                    truth_howitt_cap3, bars_howitt_cap3, weights_howitt_cap3, 'howitt_capture3', 'beacon'),
  score_assignments('results/beacon/howitt_batch1_c1/assignments.csv',
                    truth_howitt_b1c1, bars_howitt_b1c1, weights_howitt_b1c1, 'howitt_batch1_c1', 'beacon'),
  score_assignments('results/beacon/howitt_batch1_c2/assignments.csv',
                    truth_howitt_b1c2, bars_howitt_b1c2, weights_howitt_b1c2, 'howitt_batch1_c2', 'beacon'),
  score_assignments('results/beacon/howitt_batch2_c1/assignments.csv',
                    truth_howitt_b2c1, bars_howitt_b2c1, weights_howitt_b2c1, 'howitt_batch2_c1', 'beacon'),
  score_assignments('results/beacon/howitt_batch2_c2/assignments.csv',
                    truth_howitt_b2c2, bars_howitt_b2c2, weights_howitt_b2c2, 'howitt_batch2_c2', 'beacon'),
  score_assignments('results/beacon/howitt_batch3_c1/assignments.csv',
                    truth_howitt_b3c1, bars_howitt_b3c1, weights_howitt_b3c1, 'howitt_batch3_c1', 'beacon'),
  score_assignments('results/beacon/howitt_batch3_c2/assignments.csv',
                    truth_howitt_b3c2, bars_howitt_b3c2, weights_howitt_b3c2, 'howitt_batch3_c2', 'beacon')
)

print(results)
write.csv(results, 'analysis/beacon_analysis/beacon_scores.csv', row.names = FALSE)
message('Done :)')
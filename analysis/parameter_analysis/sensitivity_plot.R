#Script to plot scores

# TO RUN:
#        conda activate demux-r
#        Rscript analysis/parameter_analysis/sensitivity_plot.R

library(ggplot2)

#load results
dm2_results_df <- read.csv('analysis/parameter_analysis/dm2_mcginnis_ms.csv')
hashsolo_df <- read.csv('analysis/parameter_analysis/hashsolo_mcginnis_ms.csv')

dir.create('analysis/parameter_analysis/figures', recursive = TRUE, showWarnings = FALSE)

metrics <- c('f1_micro', 'f1_macro', 'f1_weighted', 'mcc', 'concordance')

# demultiplex2 heatmaps
for (metric in metrics) {
  p <- ggplot(dm2_results_df, aes(x = factor(init_cos_cut), y = factor(prob_cut), fill = .data[[metric]])) +
    geom_tile() +
    geom_text(aes(label = round(.data[[metric]], 3)), size = 3) +
    scale_fill_gradient2(low = 'red', mid = 'yellow', high = 'green',
                         midpoint = 0.5, limits = c(0, 1)) +
    labs(title = paste('deMULTIplex2 -', metric, '- mcginnis_ms'),
         x = 'init_cos_cut', y = 'prob_cut') +
    theme_minimal()
  
  ggsave(paste0('analysis/parameter_analysis/figures/dm2_', metric, '_mcginnis_ms.pdf'), p)
  print(p)
}

# hashsolo heatmaps - facet by n_noise
for (metric in metrics) {
  p <- ggplot(hashsolo_df, aes(x = factor(prior_neg), y = factor(prior_singlet), fill = .data[[metric]])) +
    geom_tile() +
    geom_text(aes(label = round(.data[[metric]], 3)), size = 3) +
    scale_fill_gradient2(low = 'red', mid = 'yellow', high = 'green',
                         midpoint = 0.5, limits = c(0, 1)) +
    facet_wrap(~ n_noise, labeller = label_both) +
    labs(title = paste('HashSolo -', metric, '- mcginnis_ms'),
         x = 'prior_neg', y = 'prior_singlet') +
    theme_minimal()
  
  ggsave(paste0('analysis/parameter_analysis/figures/hashsolo_', metric, '_mcginnis_ms.pdf'), p)
  print(p)
}

message(paste('Done :)'))
#Script to plot scores

# TO RUN:
#        conda activate demux-r
#        Rscript analysis/parameter_analysis/sensitivity_plot.R

library(ggplot2)
library(tidyverse)

#load results
dm2_results_df <- read.csv('analysis/parameter_analysis/dm2_mcginnis_ms.csv')
hashsolo_df <- read.csv('analysis/parameter_analysis/hashsolo_mcginnis_ms.csv')

dir.create('analysis/parameter_analysis/figures', recursive = TRUE, showWarnings = FALSE)

metrics <- c('f1_micro', 'f1_macro', 'f1_weighted', 'mcc', 'concordance')

#dm2 heatmaps
for (metric in metrics) {
  p <- ggplot(dm2_results_df, aes(x = factor(init_cos_cut), y = factor(prob_cut), fill = .data[[metric]])) +
    geom_tile() +
    geom_text(aes(label = round(.data[[metric]], 3)), size = 3) +
    scale_fill_gradient2(low = 'red3', mid = 'gold3', high = 'springgreen3',
                         midpoint = 0.5, limits = c(0, 1)) +
    labs(title = paste('deMULTIplex2 -', metric, '- mcginnis_ms'),
         x = 'init_cos_cut', y = 'prob_cut') +
    theme_minimal()
  
  ggsave(paste0('analysis/parameter_analysis/figures/dm2_', metric, '_mcginnis_ms.pdf'), p)
  print(p)
}

# dm2 line plots
dm2_long <- dm2_results_df %>%
  pivot_longer(cols = c(f1_micro, f1_macro, f1_weighted, mcc, concordance),
               names_to = 'metric', values_to = 'score')

p <- ggplot(dm2_long, aes(x = prob_cut, y = score, color = metric, group = metric)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ init_cos_cut, labeller = label_both) +
  labs(title = 'deMULTIplex2 - all metrics vs prob_cut by init_cos_cut',
       x = 'prob_cut', y = 'score') +
  theme_minimal()

ggsave('analysis/parameter_analysis/figures/dm2_lineplots_mcginnis_ms.pdf', p, width = 12, height = 8)
print(p)

#hsolo heatmaps
for (metric in metrics) {
  p <- ggplot(hashsolo_df, aes(x = factor(prior_neg), y = factor(prior_singlet), fill = .data[[metric]])) +
    geom_tile() +
    geom_text(aes(label = round(.data[[metric]], 3)), size = 3) +
    scale_fill_gradient2(low = 'red3', mid = 'gold3', high = 'springgreen3',
                         midpoint = 0.5, limits = c(0, 1)) +
    facet_wrap(~ n_noise, labeller = label_both) +
    labs(title = paste('HashSolo -', metric, '- mcginnis_ms'),
         x = 'prior_neg', y = 'prior_singlet') +
    theme_minimal()
  
  ggsave(paste0('analysis/parameter_analysis/figures/hashsolo_', metric, '_mcginnis_ms.pdf'), p)
  print(p)
}

#hsolo line plots
hashsolo_long <- hashsolo_df %>%
  pivot_longer(cols = c(f1_micro, f1_macro, f1_weighted, mcc, concordance),
               names_to = 'metric', values_to = 'score')

p <- ggplot(hashsolo_long, aes(x = prior_neg, y = score, color = metric, group = metric)) +
  geom_line() +
  geom_point() +
  facet_grid(n_noise ~ prior_singlet, labeller = label_both) +
  labs(title = 'HashSolo - all metrics vs prior_neg',
       x = 'prior_neg', y = 'score') +
  theme_minimal()

ggsave('analysis/parameter_analysis/figures/hashsolo_lineplots_mcginnis_ms.pdf', p, width = 14, height = 10)
print(p)

message(paste('Done :)'))
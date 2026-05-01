#script to run deMULTIplex

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/demultiplex/run.R dataset_# data/file_name.csv

#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

library(deMULTIplex)
library(tidyverse)

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')] #will need to check other datasets for diff col names !!

data <- as.matrix(data)

#quantile sweep across thresholds
data.sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by = 0.02)) {
  n <- n + 1
  data.sweep.list[[n]] <- classifyCells(data, q = q)
  names(data.sweep.list)[n] <- paste0("q=", q)
}

#find best quantile
best.q <- 0.51
print(paste("Using quantile:", best.q))

#get classifications
final.calls <- classifyCells(data, q = best.q)

classifications <- data.frame(
  cell_barcode = names(final.calls),
  classification = case_when(
    final.calls == 'Doublet' ~ 'multiplet',
    final.calls == 'Negative' ~ 'negative',
    TRUE ~ 'singlet'
  )
)

#save classifications
dir.create(paste0('results/demultiplex/', dataset_id), recursive = TRUE, showWarnings = FALSE)
write.csv(classifications,
          paste0('results/demultiplex/', dataset_id, '/classifications.csv'),
          row.names = FALSE)

#save summary counts
summary_counts <- classifications %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'demultiplex')

totals <- summary_counts %>%
  summarise(classification = "total", n = sum(n), dataset = dataset_id, method = 'demultiplex')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts, 
          paste0('results/demultiplex/', dataset_id, '/summary.csv'), 
          row.names = FALSE)

#move plots to results folder
if (file.exists('Rplots.pdf')) {
  file.rename('Rplots.pdf', paste0('results/demultiplex/', dataset_id, "/Rplots.pdf"))
}

#print classification counts
print(summary_counts)

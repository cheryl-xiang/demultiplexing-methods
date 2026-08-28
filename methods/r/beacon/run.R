# script to run Beacon

# to run in terminal:
#    (1) conda activate demux-r
#    (2) Rscript methods/r/beacon/run.R dataset data/dataset/hto/file_name.csv [switch_transpose]
#    switch_transpose: TRUE to switch default transposing behavior (e.g. gaublomme where barcodes are cols)

source('methods/r/beacon/beacon_suite.R')
library(tidyverse)

#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

if (length(args) >= 3) {
  switch_transpose <- as.logical(args[3])
} else {
  switch_transpose <- FALSE
}

#assign params
params <- list()
params$sig_max <- 3
params$min_count <- 2
params$max_bf <- 0.97
params$max_iter <- 30

# data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total', 'TSNE1', 'TSNE2')]

if (switch_transpose) {
  data <- as.matrix(data)
} else {
  data <- t(as.matrix(data))
}

print(paste('HTOs:', nrow(data)))
print(paste('Cells:', ncol(data)))

#run beacon
calls <- beacon_calls(data, params)

#save assignments
calls_df <- data.frame(
  cell_barcode = names(calls),
  assignment = as.numeric(calls)
)

dir.create(paste0('results/beacon/', dataset_id), recursive = TRUE, showWarnings = FALSE)
write.csv(calls_df, paste0('results/beacon/', dataset_id, '/assignments.csv'), row.names = FALSE)

#save summary counts
summary_counts <- calls_df %>%
  mutate(classification = case_when(
    assignment == 0 ~ 'negative',
    assignment == 1000 ~ 'doublet',
    TRUE ~ 'singlet'
  )) %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'beacon')

totals <- summary_counts %>%
  summarise(classification = 'total', n = sum(n), dataset = dataset_id, method = 'beacon')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts, paste0('results/beacon/', dataset_id, '/summary.csv'), row.names = FALSE)
print(summary_counts)
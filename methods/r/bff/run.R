# script to run BFF_Raw and BFF_Cluster

# to run in terminal:
#    (1) conda activate demux-r
#    (2) Rscript methods/r/bff/run.R dataset data/dataset/hto/file_name.csv

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

library(cellhashR)
library(tidyverse)

# data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')]
mat <- t(data)
mat <- mat[rowSums(mat) > 0, ]

# run BFF_Raw and BFF_Cluster together
res <- tryCatch({
  GenerateCellHashingCalls(
    mat,
    methods = c('bff_raw', 'bff_cluster'),
  )
}, error = function(e) {
  message(paste('BFF failed:', e$message))
  dir.create(paste0('results/bffraw/', dataset_id), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0('results/bffcluster/', dataset_id), recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(error = e$message),
            paste0('results/bffraw/', dataset_id, '/error.csv'), row.names = FALSE)
  write.csv(data.frame(error = e$message),
            paste0('results/bffcluster/', dataset_id, '/error.csv'), row.names = FALSE)
  NULL
})

if (is.null(res)) {
  quit(status = 0)
}

# --- BFF_Raw ---
classifications_raw <- data.frame(
  cell_barcode = res$cellbarcode,
  classification = case_when(
    res$bff_raw == 'Doublet' ~ 'multiplet',
    res$bff_raw == 'Negative' ~ 'negative',
    TRUE ~ 'singlet'
  )
)

dir.create(paste0('results/bffraw/', dataset_id), recursive = TRUE, showWarnings = FALSE)

write.csv(classifications_raw,
          paste0('results/bffraw/', dataset_id, '/classifications.csv'),
          row.names = FALSE)

summary_raw <- classifications_raw %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'bff_raw')

totals_raw <- summary_raw %>%
  summarise(classification = 'total', n = sum(n), dataset = dataset_id, method = 'bff_raw')

summary_raw <- bind_rows(summary_raw, totals_raw)

write.csv(summary_raw,
          paste0('results/bffraw/', dataset_id, '/summary.csv'),
          row.names = FALSE)

# --- BFF_Cluster ---
classifications_cluster <- data.frame(
  cell_barcode = res$cellbarcode,
  classification = case_when(
    res$bff_cluster == 'Doublet' ~ 'multiplet',
    res$bff_cluster == 'Negative' ~ 'negative',
    TRUE ~ 'singlet'
  )
)

dir.create(paste0('results/bffcluster/', dataset_id), recursive = TRUE, showWarnings = FALSE)

write.csv(classifications_cluster,
          paste0('results/bffcluster/', dataset_id, '/classifications.csv'),
          row.names = FALSE)

summary_cluster <- classifications_cluster %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'bff_cluster')

totals_cluster <- summary_cluster %>%
  summarise(classification = 'total', n = sum(n), dataset = dataset_id, method = 'bff_cluster')

summary_cluster <- bind_rows(summary_cluster, totals_cluster)

write.csv(summary_cluster,
          paste0('results/bffcluster/', dataset_id, '/summary.csv'),
          row.names = FALSE)

# move plots
if (file.exists('Rplots.pdf')) {
  file.rename('Rplots.pdf', paste0('results/bffraw/', dataset_id, '/Rplots.pdf'))
}

# print results
print("BFF_Raw:")
print(summary_raw)
print("BFF_Cluster:")
print(summary_cluster)
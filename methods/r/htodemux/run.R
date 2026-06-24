#script to run HTODemux

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/htodemux/run.R dataset data/dataset/hto/file_name.csv [n_barcodes] [switch_transpose]
#    switch_transpose: TRUE to switch default transposing behavior (where barcodes are cols)

args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

if (length(args) >= 3) {
  if (args[3] %in% c('TRUE', 'FALSE', 'true', 'false')) {
    n_barcodes <- NULL
    switch_transpose <- as.logical(args[3])
  } else {
    n_barcodes <- as.integer(args[3])
    if (length(args) >= 4) {
      switch_transpose <- as.logical(args[4])
    } else {
      switch_transpose <- FALSE
    }
  }
} else {
  n_barcodes <- NULL
  switch_transpose <- FALSE
}

library(Seurat)
library(tidyverse)

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')]

#select n_barcodes if specified
if (!is.null(n_barcodes)) {
  data <- data[, 1:n_barcodes]
}

#default is to transpose (cells as rows -> barcodes as rows, cells as cols)
#switch_transpose = TRUE flips this for datasets like gaublomme
do_transpose <- !switch_transpose

if (do_transpose) {
  mat <- t(as.matrix(data))
} else {
  mat <- as.matrix(data)
}

print(paste('Barcodes:', nrow(mat)))
print(paste('Cells:', ncol(mat)))

mat_sparse <- Matrix::Matrix(mat, sparse = TRUE)

seurat_obj <- CreateSeuratObject(counts = mat_sparse)
seurat_obj[['HTO']] <- CreateAssayObject(counts = mat_sparse)
seurat_obj <- NormalizeData(seurat_obj, assay = 'HTO', normalization.method = 'CLR')
res <- HTODemux(seurat_obj, assay = 'HTO', positive.quantile = 0.99, kfunc = 'kmeans')

classifications <- data.frame(
  cell_barcode = colnames(seurat_obj),
  classification = case_when(
    res$HTO_classification.global == 'Doublet' ~ 'doublet',
    res$HTO_classification.global == 'Negative' ~ 'negative',
    TRUE ~ 'singlet'
  )
)

dir.create(paste0('results/htodemux/', dataset_id), recursive = TRUE, showWarnings = FALSE)
write.csv(classifications,
          paste0('results/htodemux/', dataset_id, '/classifications.csv'),
          row.names = FALSE)

summary_counts <- classifications %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'htodemux')

totals <- summary_counts %>%
  summarise(classification = 'total', n = sum(n), dataset = dataset_id, method = 'htodemux')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts,
          paste0('results/htodemux/', dataset_id, '/summary.csv'),
          row.names = FALSE)

print(summary_counts)
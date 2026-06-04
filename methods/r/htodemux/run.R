#script to run HTODemux

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/htodemux/run.R dataset data/dataset/hto/file_name.csv n_barcodes

#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]
if (length(args) >= 3) {
  n_barcodes <- as.integer(args[3])
} else {
  n_barcodes <- NULL
}

library(Seurat)
library(tidyverse) 

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')]

if (!is.null(n_barcodes)) {
  data <- data[, 1:n_barcodes]
  print(paste('Using first', n_barcodes, 'barcodes'))
} else {
  print(paste('Using all', ncol(data), 'barcodes'))
}

#transpose
mat <- t(as.matrix(data))
mat_sparse <- Matrix::Matrix(mat, sparse = TRUE)

seurat_obj <- CreateSeuratObject(counts = mat_sparse)

#add barcode data as a new assay independent from RNA
seurat_obj[["HTO"]] <- CreateAssayObject(counts = mat_sparse)

#normalize
seurat_obj <- NormalizeData(seurat_obj, assay = 'HTO', normalization.method = "CLR")

#run HTOdemux
res <-  HTODemux(seurat_obj, assay = 'HTO', positive.quantile = 0.99)
#res <-  HTODemux(seurat_obj, assay = 'HTO', positive.quantile = 0.99, kfunc = 'kmeans')  #<-- Runs on Winkler, but not default


#get classifications
classifications <- data.frame(
  cell_barcode = colnames(seurat_obj),
  classification = case_when(
    res$HTO_classification.global == 'Doublet' ~ 'doublet',
    res$HTO_classification.global == 'Negative' ~ 'negative',
    TRUE ~ 'singlet'
  )
)

#save classifications
dir.create(paste0('results/htodemux/', dataset_id), recursive = TRUE, showWarnings = FALSE)
write.csv(classifications,
          paste0('results/htodemux/', dataset_id, '/classifications.csv'),
          row.names = FALSE)

#save summary counts
summary_counts <- classifications %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'htodemux')

totals <- summary_counts %>%
  summarise(classification = "total", n = sum(n), dataset = dataset_id, method = 'htodemux')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts, 
          paste0('results/htodemux/', dataset_id, '/summary.csv'), 
          row.names = FALSE)

#print classification counts
print(summary_counts)

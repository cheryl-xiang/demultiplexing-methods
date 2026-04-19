#script to run HTODemux

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/htodemux/run.R dataset_# data/file_name.csv

#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

library(Seurat)
library(tidyverse) 

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')]  #will need to check other datasets for diff col names !!

mat <- t(data)     #expects barcodes as rows, cells as cols

seurat_obj <- CreateSeuratObject(counts = mat)

#add barcode data as a new assay independent from RNA
seurat_obj[['barcode']] <- CreateAssayObject(counts = mat)

#normalize
seurat_obj <- NormalizeData(seurat_obj, assay = 'barcode', normalization.method = "CLR")

#run HTOdemux
res <-  HTODemux(seurat_obj, assay = 'barcode', positive.quantile = 0.99)

#get classifications
classifications <- data.frame(
  cell_barcode = colnames(seurat_obj),
  classification = case_when(
    res$barcode_classification.global == "Doublet" ~ "multiplet",
    res$barcode_classification.global == "Negative" ~ "negative",
    TRUE ~ "singlet"
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



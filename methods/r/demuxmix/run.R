#script to run DemuxMix

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/demuxmix/run.R dataset_# data/file_name.csv

#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

library(demuxmix)
library(tidyverse)

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')] #will need to check other datasets for diff col names !!

mat <- t(data) #expects cells as cols, barcodes as rows
cell_totals <- colSums(mat)

#remove empty droplets
mat <- mat[, cell_totals >= 1]
cell_totals <- cell_totals[cell_totals >= 1]

#run demuxmix- 
res <- demuxmix(mat, rna = cell_totals)

classes <- dmmClassify(res)

#get classifications
classifications <- data.frame(
  cell_barcode = rownames(classes),
  classification = case_when(
    classes$Type == 'multiplet' ~ 'multiplet',
    classes$Type == 'singlet' ~ 'singlet',
    TRUE ~ 'negative'
  )
)

#save classifications
dir.create(paste0('results/demuxmix/', dataset_id), recursive = TRUE, showWarnings = FALSE)
write.csv(classifications,
          paste0('results/demuxmix/', dataset_id, '/classifications.csv'),
          row.names = FALSE)

#save summary counts
summary_counts <- classifications %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'demuxmix')

totals <- summary_counts %>%
  summarise(classification = "total", n = sum(n), dataset = dataset_id, method = 'demuxmix')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts, 
          paste0('results/demuxmix/', dataset_id, '/summary.csv'), 
          row.names = FALSE)

#move plots to results folder
if (file.exists('Rplots.pdf')) {
  file.rename('Rplots.pdf', paste0('results/demuxmix/', dataset_id, "/Rplots.pdf"))
}

#print classification counts
print(summary_counts)
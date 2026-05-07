#script to run DemuxMix Naive

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/demuxmixnaive/run.R dataset_# data/dataset_#/hto/file_name.csv data/dataset_#/rna/

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

#run demuxmix naive
res <- demuxmix(mat, model = 'naive')
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
dir.create(paste0('results/demuxmixnaive/', dataset_id), recursive = TRUE, showWarnings = FALSE)
write.csv(classifications,
          paste0('results/demuxmixnaive/', dataset_id, '/classifications.csv'),
          row.names = FALSE)

#save summary counts
summary_counts <- classifications %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'demuxmix naive')

totals <- summary_counts %>%
  summarise(classification = "total", n = sum(n), dataset = dataset_id, method = 'demuxmix naive')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts, 
          paste0('results/demuxmixnaive/', dataset_id, '/summary.csv'), 
          row.names = FALSE)

#move plots to results folder
if (file.exists('Rplots.pdf')) {
  file.rename('Rplots.pdf', paste0('results/demuxmixnaive/', dataset_id, "/Rplots.pdf"))
}

#print classification counts
print(summary_counts)
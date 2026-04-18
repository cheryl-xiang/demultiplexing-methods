#script to run HashedDrops

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/hasheddrops/run.R dataset_# data/file_name.csv

#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

library(DropletUtils)
library(tidyverse)   #hmm check if you need this

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')]  #will need to check other datasets for diff col names !!

mat <- t(data)     #expects barcodes as rows, cells as cols

#run HashedDrops
res <- hashedDrops(mat)

#get classifications
classifications <- data.frame(
  cell_barcode = colnames(mat),
  classification = case_when(
    res$Doublet == TRUE ~ 'multiplet',
    res$Confident == TRUE ~ "singlet",
    TRUE ~ "negative"
  )
)

#save classifications
dir.create(paste0('results/hasheddrops/', dataset_id), recursive = TRUE, showWarnings = FALSE)
write.csv(classifications,
          paste0('results/hasheddrops/', dataset_id, '/classifications.csv'),
          row.names = FALSE)

#save summary counts
summary_counts <- classifications %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'hasheddrops')

totals <- summary_counts %>%
  summarise(classification = "total", n = sum(n), dataset = dataset_id, method = 'hasheddrops')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts, 
          paste0('results/hasheddrops/', dataset_id, '/summary.csv'), 
          row.names = FALSE)

#move plots to results folder
if (file.exists('Rplots.pdf')) {
  file.rename('Rplots.pdf', paste0('results/hasheddrops/', dataset_id, "/Rplots.pdf"))
}

#print classification counts
print(summary_counts)




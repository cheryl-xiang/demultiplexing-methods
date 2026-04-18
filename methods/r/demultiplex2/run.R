#script to run deMULTIplex2

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/demultiplex2/run.R dataset_# data/file_name.csv

#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

library(deMULTIplex2)
library(tidyverse)

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')] #will need to check other datasets for diff col names !!

#run deMULTIplex2
res <- demultiplexTags(data)

#get classifications
classifications <- res$assign_table %>%
  rownames_to_column('cell_barcode') %>%
  select(cell_barcode, classification = final_assign)   #ask if i should just save as multiplet/singlet/negative rather than specifc barcode

#save classifications
dir.create(paste0('results/demultiplex2/', dataset_id), recursive = TRUE, showWarnings = FALSE)
write.csv(classifications,
          paste0('results/demultiplex2/', dataset_id, '/classifications.csv'),
          row.names = FALSE)

#save summary counts
summary_counts <- classifications %>%
  mutate(classification = case_when(
    classification == 'negative' ~ 'negative',
    str_detect(classification, 'multiplet') ~ 'multiplet',
    TRUE ~ 'singlet'
  )) %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'demultiplex2')

totals <- summary_counts %>%
  summarise(classification = 'total', n = sum(n), dataset = dataset_id, method = 'demultiplex2')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts, 
          paste0('results/demultiplex2/', dataset_id, '/summary.csv'), 
          row.names = FALSE)

#move assignment pdf to results folder
pdf_file <- list.files(pattern = ".*assignment\\.pdf$")
if (length(pdf_file) > 0) {
  file.rename(pdf_file, paste0('results/demultiplex2/', dataset_id, '/', pdf_file))
}

#move Rplots.pdf to results folder
if (file.exists('Rplots.pdf')) {
 file.rename('Rplots.pdf', paste0('results/demultiplex2/', dataset_id, '/Rplots.pdf'))
}

#print classification counts
print(summary_counts)


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
data <- read.csv('data/GSM4904942_8donor_PBMC_AH_MULTI_matrix.csv', row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')]

#run deMULTIplex2
res <- demultiplexTags(data)

#get classifications
classifications <- res$assign_table %>%
  rownames_to_column('cell_barcode') %>%
  select(cell_barcode, classification = final_assign)

#save classifications
dir.create('results/demultiplex2/dataset_1', recursive = TRUE, showWarnings = FALSE)
write.csv(classifications,
          'results/demultiplex2/dataset_1/classifications.csv',
          row.names = FALSE)

#save summary counts
summary_counts <- classifications %>%
  mutate(classification = case_when(
    classification == 'negative' ~ 'negative',
    str_detect(classification, 'multiplet') ~ 'multiplet',
    TRUE ~ 'singlet'
  )) %>%
  count(classification) %>%
  mutate(dataset = 'dataset_1', method = 'demultiplex2')

write.csv(summary_counts,
          'results/demultiplex2/dataset_1/summary.csv',
          row.names = FALSE)

# move assignment pdf to results folder
pdf_file <- list.files(pattern = ".*assignment\\.pdf$")
if (length(pdf_file) > 0) {
  file.rename(pdf_file, paste0('results/demultiplex2/dataset_1/', pdf_file))
}

# move Rplots.pdf to results folder
if (file.exists('Rplots.pdf')) {
  file.rename('Rplots.pdf', 'results/demultiplex2/dataset_1/Rplots.pdf')
}

#print classification counts
print(summary_counts)


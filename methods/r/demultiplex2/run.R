#script to run deMULTIplex2

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/demultiplex2/run.R dataset data/dataset/hto/file_name.csv [switch_transpose]
#    switch_transpose: TRUE to switch default transposing behavior (where barcodes are cols)


#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

if (length(args) >= 3) {
  switch_transpose <- as.logical(args[3])
} else {
  switch_transpose <- FALSE
}

library(deMULTIplex2)
library(tidyverse)

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total', 'TSNE1', 'TSNE2')]

if (switch_transpose) {
  data <- t(data)
}

#run deMULTIplex2
res <- demultiplexTags(data)

#get classifications
classifications <- res$assign_table %>%
  rownames_to_column('cell_barcode') %>%
  select(cell_barcode, classification = final_assign) 

#add dropped cells as negative
all_cells <- rownames(data)
missing_cells <- setdiff(all_cells, classifications$cell_barcode)
if (length(missing_cells) > 0) {
  missing_df <- data.frame(
    cell_barcode = missing_cells,
    classification = "negative"
  )
  classifications <- bind_rows(classifications, missing_df)
}



#save classifications
dir.create(paste0('results/demultiplex2/', dataset_id), recursive = TRUE, showWarnings = FALSE)
write.csv(classifications,
          paste0('results/demultiplex2/', dataset_id, '/classifications.csv'),
          row.names = FALSE)

#save summary counts
summary_counts <- classifications %>%
  mutate(classification = case_when(
    classification == 'negative' ~ 'negative',
    str_detect(classification, 'multiplet') ~ 'doublet',
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


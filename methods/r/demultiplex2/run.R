#script to run deMULTIplex2

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/demultiplex2/run.R dataset data/dataset/hto/file_name.csv [switch_transpose] [barcode_map]
#    switch_transpose: TRUE to switch default transposing behavior (e.g. gaublomme)
#    barcode_map: optional csv with columns 'barcode' and 'index' for numeric assignment

args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

if (length(args) >= 3) {
  switch_transpose <- as.logical(args[3])
} else {
  switch_transpose <- FALSE
}

barcode_map <- NULL
if (length(args) >= 4) {
  barcode_map <- read.csv(args[4])
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

#create output directory
dir.create(paste0('results/demultiplex2/', dataset_id), recursive = TRUE, showWarnings = FALSE)

#get raw assignments
raw_assignments <- res$assign_table %>%
  rownames_to_column('cell_barcode') %>%
  select(cell_barcode, final_assign)

#map to numeric assignments
if (!is.null(barcode_map)) {
  barcode_lookup <- setNames(barcode_map$index, barcode_map$barcode)
} else {
  barcode_lookup <- c()
}

assignments <- raw_assignments %>%
  mutate(assignment = case_when(
    str_detect(final_assign, 'negative') ~ 0,
    str_detect(final_assign, 'multiplet') ~ 1000,
    final_assign %in% names(barcode_lookup) ~ barcode_lookup[final_assign],
    TRUE ~ NA_real_
  )) %>%
  select(cell_barcode, assignment)

write.csv(assignments,
          paste0('results/demultiplex2/', dataset_id, '/assignments.csv'),
          row.names = FALSE)

#save summary counts
summary_counts <- assignments %>%
  mutate(classification = case_when(
    assignment == 0 ~ 'negative',
    assignment == 1000 ~ 'doublet',
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

#move pdfs
pdf_file <- list.files(pattern = '.*assignment\\.pdf$')
if (length(pdf_file) > 0) {
  file.rename(pdf_file, paste0('results/demultiplex2/', dataset_id, '/', pdf_file))
}

if (file.exists('Rplots.pdf')) {
  file.rename('Rplots.pdf', paste0('results/demultiplex2/', dataset_id, '/Rplots.pdf'))
}

print(summary_counts)
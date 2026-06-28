#script to run HashedDrops

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/hasheddrops/run.R dataset data/dataset/hto/file_name.csv [switch_transpose]
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

library(DropletUtils)
library(tidyverse)   #hmm check if you need this

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total', 'TSNE1', 'TSNE2')]

do_transpose <- !switch_transpose

if (do_transpose) {
  mat <- t(as.matrix(data))
} else {
  mat <- as.matrix(data)
}

#run HashedDrops
res <- hashedDrops(mat)

#get classifications
classifications <- data.frame(
  cell_barcode = colnames(mat),
  classification = case_when(
    res$Doublet == TRUE ~ 'doublet',
    res$Confident == TRUE ~ "singlet",
    TRUE ~ "negative"
  )
)

#barcode calls
barcode_calls <- data.frame(
  cell_barcode = colnames(mat),
  call = res$Best  #barcode index for all cells by default
)
barcode_calls$call[res$Confident == FALSE] <- 0    # negative
barcode_calls$call[res$Doublet == TRUE] <- 1000    # doublet

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




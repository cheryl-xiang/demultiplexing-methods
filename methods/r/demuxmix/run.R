#script to run DemuxMix

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/demuxmix/run.R dataset data/dataset/hto/file_name.csv data/dataset/rna/file.rds [switch_transpose]
#    switch_transpose: TRUE to switch default transposing behavior (where barcodes are cols)


#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]
rna_dir <- args[3]

if (length(args) >= 4) {
  switch_transpose <- as.logical(args[4])
} else {
  switch_transpose <- FALSE
}

library(demuxmix)
library(tidyverse)
library(Matrix)

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total', 'TSNE1', 'TSNE2')]

do_transpose <- !switch_transpose

if (do_transpose) {
  mat <- t(as.matrix(data))
} else {
  mat <- as.matrix(data)
} #expects cells as cols, barcodes as rows

#load rna data
if (grepl('\\.rds$', rna_dir, ignore.case = TRUE)) {
  rna_mat <- readRDS(rna_dir)
} else if (grepl('\\.h5$', rna_dir, ignore.case = TRUE)) {
  library(Seurat)
  rna_mat <- Read10X_h5(rna_dir)
} else {
  barcodes <- read.table(file.path(rna_dir, 'barcodes.tsv'), header = FALSE)$V1
  features <- read.table(file.path(rna_dir, 'features.tsv'), header = FALSE)
  rna_mat <- readMM(file.path(rna_dir, 'matrix.mtx'))
  rownames(rna_mat) <- features$V2
  colnames(rna_mat) <- barcodes
}

#strip -1 suffix from RNA barcodes
colnames(rna_mat) <- sub('-1$', '', colnames(rna_mat))

#find common cells
common_cells <- intersect(colnames(mat), colnames(rna_mat))
print(paste('HTO cells:', ncol(mat)))
print(paste('RNA cells:', ncol(rna_mat)))
print(paste('Common cells:', length(common_cells)))
mat <- mat[, common_cells, drop = FALSE]
mat <- mat[rowSums(mat) > 0, , drop = FALSE]
rna_counts <- Matrix::colSums(rna_mat[, common_cells] > 0)

#run demuxmix
res <- demuxmix(mat, rna = rna_counts)
classes <- dmmClassify(res)

#get classifications
classifications <- data.frame(
  cell_barcode = rownames(classes),
  classification = case_when(
    classes$Type == 'multiplet' ~ 'doublet',
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
  summarise(classification = 'total', n = sum(n), dataset = dataset_id, method = 'demuxmix')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts, 
          paste0('results/demuxmix/', dataset_id, '/summary.csv'), 
          row.names = FALSE)

#move plots to results folder
if (file.exists('Rplots.pdf')) {
  file.rename('Rplots.pdf', paste0('results/demuxmix/', dataset_id, '/Rplots.pdf'))
}

#print classification counts
print(summary_counts)
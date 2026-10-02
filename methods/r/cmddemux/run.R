# script to run CMDdemux

# to run in terminal:
#    (1) conda activate demux-cmd
#    (2) Rscript methods/r/cmddemux/run.R dataset data/dataset/hto/file_name.csv [switch_transpose] [barcode_map]
#    switch_transpose: TRUE to switch default transposing behavior (e.g. gaublomme)
#    barcode_map: optional csv with columns 'barcode' and 'index' for numeric assignment

#IQR function was written wrong in the package? 
iqr <- function(x) IQR(x)

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

library(CMDdemux)
library(tidyverse)

# data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total', 'TSNE1', 'TSNE2')]

# CMDdemux expects barcodes as rows, cells as cols
if (switch_transpose) {
  hash.count <- as.matrix(data)
} else {
  hash.count <- t(as.matrix(data))
}

print(paste('HTOs:', nrow(hash.count)))
print(paste('Cells:', ncol(hash.count)))

# create output directory
dir.create(paste0('results/cmddemux/', dataset_id), recursive = TRUE, showWarnings = FALSE)

# run CMDdemux with timing
start_time <- proc.time()

clr.norm <- LocalCLRNorm(hash.count)
kmed.cl <- KmedCluster(clr.norm)
cl.dist <- EuclideanClusterDist(clr.norm, kmed.cl)
noncore <- DefineNonCore(cl.dist, kmed.cl)
cluster.assign <- LabelClusterHTO(clr.norm, kmed.cl, noncore, label_method = 'medoids')
md.mat <- CalculateMD(clr.norm, noncore, kmed.cl, cluster.assign)
outlier.assign <- AssignOutlierDrop(md.mat)
cmddemux.assign <- CMDdemuxClass(md.mat, hash.count, outlier.assign, use_gex_data = FALSE)

elapsed <- proc.time() - start_time
runtime <- as.numeric(elapsed['elapsed'])

#save runtime
runtime_df <- data.frame(
  dataset = dataset_id,
  method = 'cmddemux',
  runtime_seconds = runtime
)
write.csv(runtime_df, paste0('results/cmddemux/', dataset_id, '/runtime.csv'), row.names = FALSE)
message(paste('Runtime:', round(runtime, 2), 'seconds'))

# map to numeric assignments
if (!is.null(barcode_map)) {
  barcode_lookup <- setNames(barcode_map$index, barcode_map$barcode)
} else {
  barcode_lookup <- c()
}

assignments <- data.frame(
  cell_barcode = colnames(md.mat),
  global_class = cmddemux.assign$global_class,
  demux_id = cmddemux.assign$demux_id
) %>%
  mutate(assignment = case_when(
    global_class == 'Negative' ~ 0,
    global_class == 'Doublet' ~ 1000,
    demux_id %in% names(barcode_lookup) ~ barcode_lookup[demux_id],
    TRUE ~ NA_real_
  )) %>%
  select(cell_barcode, assignment)

write.csv(assignments,
          paste0('results/cmddemux/', dataset_id, '/assignments.csv'),
          row.names = FALSE)

# save summary counts
summary_counts <- assignments %>%
  mutate(classification = case_when(
    assignment == 0 ~ 'negative',
    assignment == 1000 ~ 'doublet',
    TRUE ~ 'singlet'
  )) %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'cmddemux')

totals <- summary_counts %>%
  summarise(classification = 'total', n = sum(n), dataset = dataset_id, method = 'cmddemux')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts,
          paste0('results/cmddemux/', dataset_id, '/summary.csv'),
          row.names = FALSE)

if (file.exists('Rplots.pdf')) {
  file.rename('Rplots.pdf', paste0('results/cmddemux/', dataset_id, '/Rplots.pdf'))
}

print(summary_counts)
#script to run deMULTIplex
#cells as rows, barcodes as columns

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/demultiplex/run.R dataset data/dataset/hto/file_name.csv n_rounds rescue_threshold [n_barcodes] [switch_transpose]
#    switch_transpose: TRUE to switch default transposing behavior (where barcodes are cols)

args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]
n_rounds <- as.integer(args[3])
rescue_threshold <- as.integer(args[4])

if (length(args) >= 5) {
  if (args[5] %in% c('TRUE', 'FALSE', 'true', 'false')) {
    n_barcodes <- NULL
    switch_transpose <- as.logical(args[5])
  } else {
    n_barcodes <- as.integer(args[5])
    if (length(args) >= 6) {
      switch_transpose <- as.logical(args[6])
    } else {
      switch_transpose <- FALSE
    }
  }
} else {
  n_barcodes <- NULL
  switch_transpose <- FALSE
}

library(deMULTIplex)
library(tidyverse)

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total', 'TSNE1', 'TSNE2')]

#select barcodes
mat <- as.matrix(data)
if (!is.null(n_barcodes)) {
  mat <- mat[, 1:n_barcodes]
  print(paste('Using first', n_barcodes, 'barcodes'))
} else {
  print(paste('Using all', ncol(mat), 'barcodes'))
}

if (switch_transpose) {
  mat <- t(mat)
}

data.full <- mat
data <- mat
neg.cells <- c()

#quantile sweeps
for (round in 1:n_rounds) {
  bar.table_sweep.list <- list()
  n <- 0
  for (q in seq(0.01, 0.99, by = 0.02)) {
    n <- n + 1
    bar.table_sweep.list[[n]] <- classifyCells(data, q = q)
    names(bar.table_sweep.list)[n] <- paste0('q=', q)
  }

  threshold.results <- tryCatch({
    findThresh(call.list = bar.table_sweep.list)
  }, error = function(e) {
    print(paste('findThresh failed in round', round, ':', e$message))
    NULL
  })

  if (is.null(threshold.results)) {
    print(paste('Stopping after round', round - 1))
    break
  }

  best.q <- findQ(threshold.results$res, threshold.results$extrema)

  if (length(best.q) == 0) {
    print(paste('No threshold found in round', round, '- stopping early'))
    break
  }

  round.calls <- classifyCells(data, q = best.q)
  new.neg.cells <- names(round.calls)[which(round.calls == 'Negative')]

  if (length(new.neg.cells) == 0) {
    print(paste('No new negatives in round', round, '- stopping early'))
    final.round.calls <- round.calls
    break
  }

  neg.cells <- unique(c(neg.cells, new.neg.cells))
  data <- data[-which(rownames(data) %in% new.neg.cells), ]
  print(paste('Round', round, 'negatives:', length(new.neg.cells)))
  print(paste('Cells remaining:', nrow(data)))

  final.round.calls <- round.calls
}

#build final calls
final.calls <- c(final.round.calls, rep('Negative', length(neg.cells)))
names(final.calls) <- c(names(final.round.calls), neg.cells)

#rescue negatives
if (rescue_threshold > 0) {
  reclass.cells <- findReclassCells(data.full, names(final.calls)[which(final.calls == 'Negative')])
  reclass.res <- rescueCells(data.full, final.calls, reclass.cells)
  rescue.ind <- which(reclass.cells$ClassStability >= rescue_threshold)
  final.calls.rescued <- final.calls
  final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]
} else {
  final.calls.rescued <- final.calls
}

#reorder rows
final.calls.rescued <- final.calls.rescued[match(rownames(data.full), names(final.calls.rescued))]

#get classifications
classifications <- data.frame(
  cell_barcode = names(final.calls.rescued),
  classification = case_when(
    final.calls.rescued == 'Doublet' ~ 'doublet',
    final.calls.rescued == 'Negative' ~ 'negative',
    TRUE ~ 'singlet'
  )
)

#save classifications
dir.create(paste0('results/demultiplex/', dataset_id), recursive = TRUE, showWarnings = FALSE)
write.csv(classifications,
          paste0('results/demultiplex/', dataset_id, '/classifications.csv'),
          row.names = FALSE)

#save summary counts
summary_counts <- classifications %>%
  count(classification) %>%
  mutate(dataset = dataset_id, method = 'demultiplex')

totals <- summary_counts %>%
  summarise(classification = 'total', n = sum(n), dataset = dataset_id, method = 'demultiplex')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts,
          paste0('results/demultiplex/', dataset_id, '/summary.csv'),
          row.names = FALSE)

if (file.exists('Rplots.pdf')) {
  file.rename('Rplots.pdf', paste0('results/demultiplex/', dataset_id, '/Rplots.pdf'))
}

print(summary_counts)
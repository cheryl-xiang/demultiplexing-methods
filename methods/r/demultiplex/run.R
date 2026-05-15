#script to run deMULTIplex

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/demultiplex/run.R dataset_# data/dataset_#/hto/file_name.csv n_rounds rescue_threshold

#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

n_rounds <- as.integer(args[3])      #number of quantile sweep rounds
rescue_threshold <- as.integer(args[4])  #class stability threshold for rescue

library(deMULTIplex)
library(tidyverse)

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')] #will need to check other datasets for diff col names !!
data.full <- as.matrix(data)  #keep full matrix for rescue step
data <- data.full

data <- as.matrix(data)

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
  threshold.results <- findThresh(call.list = bar.table_sweep.list)
  round.calls <- classifyCells(data, q = findQ(threshold.results$res, threshold.results$extrema))
  new.neg.cells <- names(round.calls)[which(round.calls == 'Negative')]
  neg.cells <- unique(c(neg.cells, new.neg.cells))
  data <- data[-which(rownames(data) %in% neg.cells), ]
  print(paste('Round', round, 'negatives:', length(new.neg.cells)))
  print(paste('Cells remaining:', nrow(data)))

  #store last round calls
  final.round.calls <- round.calls
}

#build final calls
final.calls <- c(final.round.calls, rep('Negative', length(neg.cells)))
names(final.calls) <- c(names(final.round.calls), neg.cells)

#rescue negatives
reclass.cells <- findReclassCells(data.full, names(final.calls)[which(final.calls == 'Negative')])
reclass.res <- rescueCells(data.full, final.calls, reclass.cells)

rescue.ind <- which(reclass.cells$ClassStability >= rescue_threshold)
final.calls.rescued <- final.calls
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]

#reorder rows
final.calls.rescued <- final.calls.rescued[match(rownames(data.full), names(final.calls.rescued))]

#get classifications
classifications <- data.frame(
  cell_barcode = names(final.calls.rescued),
  classification = case_when(
    final.calls.rescued == 'Doublet' ~ 'multiplet',
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
  summarise(classification = "total", n = sum(n), dataset = dataset_id, method = 'demultiplex')

summary_counts <- bind_rows(summary_counts, totals)

write.csv(summary_counts, 
          paste0('results/demultiplex/', dataset_id, '/summary.csv'), 
          row.names = FALSE)

#move plots to results folder
if (file.exists('Rplots.pdf')) {
  file.rename('Rplots.pdf', paste0('results/demultiplex/', dataset_id, "/Rplots.pdf"))
}

#print classification counts
print(summary_counts)

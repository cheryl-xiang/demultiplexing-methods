#script to run BFF_Raw

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/bffraw/run.R dataset_# data/file_name.csv

#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

library(cellhashR)
library(tidyverse)

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')] #will need to check other datasets for diff col names !!

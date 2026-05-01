#script to run deMULTIplex

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/demultiplex/run.R dataset_# data/file_name.csv

#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]

library(deMULTIplex)
library(tidyverse)

#data loading
data <- read.csv(input_file, row.names = 1)
data <- data[, !colnames(data) %in% c('nUMI', 'nUMI_total')] #will need to check other datasets for diff col names !!

#run deMULTIplex
res <- demultiplexTags(data)
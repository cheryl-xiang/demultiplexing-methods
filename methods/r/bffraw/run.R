#script to run BFF_Raw

#to run in terminal: 
#    (1) conda activate demux-r 
#    (2) Rscript methods/r/bffraw/run.R dataset_# data/file_name.csv

#read command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[1]
input_file <- args[2]
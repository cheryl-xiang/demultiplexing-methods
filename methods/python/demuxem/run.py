#script to run DemuxEM

#to run in terminal: 
#    (1) conda activate demux-py
#    (2) python3 methods/python/demuxem/run.py dataset_# data/file_name.csv

import sys
import anndata
import pandas as pd
import scanpy.external as sce
import os

#read command line arguments
dataset_id = sys.argv[1]
input_file = sys.argv[2]

#data loading
data = pd.read_csv(input_file, index_col=0)
data = data.drop(columns=[col for col in data.columns if 'nUMI' in col])  #also check for other col names in other data

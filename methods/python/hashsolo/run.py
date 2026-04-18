#script to run HashSolo

#to run in terminal: 
#    (1) conda activate demux-py
#    (2) python3 methods/python/hashsolo/run.py dataset_# data/file_name.csv

#python3 methods/python/hashsolo/run.py dataset_1 data/GSM4904942_8donor_PBMC_AH_MULTI_matrix.csv

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

adata = anndata.AnnData(X=data.values, 
                         obs=pd.DataFrame(index=data.index),
                         var=pd.DataFrame(index=data.columns))

#get hashtag column names
for col in data.columns:
    adata.obs[col] = data[col].values

hashtag_cols = list(data.columns)

#run HashSolo
sce.pp.hashsolo(adata, hashtag_cols)

#extract classifications
classifications = adata.obs[['most_likely_hypothesis']].copy()
classifications.index.name = 'cell_barcode'
classifications = classifications.reset_index()
classifications.columns = ['cell_barcode', 'classification']

#match classfication names 
classifications['classification'] = classifications['classification'].apply(
    lambda x: 'negative' if x == 0 else ('multiplet' if x == 2 else 'singlet')
)

#save classifications
os.makedirs(f'results/hashsolo/{dataset_id}', exist_ok=True)

classifications.to_csv(
    f'results/hashsolo/{dataset_id}/classifications.csv',
    index=False
)

#save summary counts
summary = classifications.groupby('classification').size().reset_index(name='n')
summary['dataset'] = dataset_id
summary['method'] = 'hashsolo'

total = pd.DataFrame([{
    'classification': 'total',
    'n': summary['n'].sum(),
    'dataset': dataset_id,
    'method': 'hashsolo'
}])

summary = pd.concat([summary, total], ignore_index=True)

summary.to_csv(
    f'results/hashsolo/{dataset_id}/summary.csv',
    index=False
)

print(summary)

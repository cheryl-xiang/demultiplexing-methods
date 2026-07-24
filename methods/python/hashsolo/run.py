#script to run HashSolo

#to run in terminal: 
#    (1) conda activate demux-py
#    (2) python3 methods/python/hashsolo/run.py dataset data/dataset/hto/file_name.csv [switch_transpose] [barcode_map]
#    switch_transpose: TRUE to switch default transposing behavior (e.g. gaublomme)
#    barcode_map: optional csv with columns 'barcode' and 'index' for numeric assignment

import sys
import anndata
import pandas as pd
import scanpy.external as sce
import os

#read command line arguments
dataset_id = sys.argv[1]
input_file = sys.argv[2]
switch_transpose = sys.argv[3].lower() == 'true' if len(sys.argv) >= 4 else False
barcode_map_file = sys.argv[4] if len(sys.argv) >= 5 else None

#data loading
data = pd.read_csv(input_file, index_col=0)
data = data.drop(columns=[col for col in data.columns if 'nUMI' in col])
data = data.drop(columns=[col for col in data.columns if 'TSNE' in col])

if switch_transpose:
    data = data.T

adata = anndata.AnnData(X=data.values,
                        obs=pd.DataFrame(index=data.index),
                        var=pd.DataFrame(index=data.columns))

for col in data.columns:
    adata.obs[col] = data[col].values

hashtag_cols = list(data.columns)

#run HashSolo
sce.pp.hashsolo(adata, hashtag_cols)

#extract classifications
classifications = adata.obs[['most_likely_hypothesis', 'Classification']].copy()
classifications.index.name = 'cell_barcode'
classifications = classifications.reset_index()

#save sample-specific numeric assignments for scoring
if barcode_map_file is not None:
    barcode_map = pd.read_csv(barcode_map_file)
    barcode_lookup = dict(zip(barcode_map['barcode'], barcode_map['index']))

    def map_to_numeric(row):
        if row['most_likely_hypothesis'] == 0:
            return 0
        elif row['most_likely_hypothesis'] == 2:
            return 1000
        else:
            return barcode_lookup.get(row['Classification'], None)

    assignments = classifications[['cell_barcode']].copy()
    assignments['assignment'] = classifications.apply(map_to_numeric, axis=1)
    assignments.to_csv(f'results/hashsolo/{dataset_id}/assignments.csv', index=False)

#standardize classifications
classifications['classification'] = classifications['most_likely_hypothesis'].apply(
    lambda x: 'negative' if x == 0 else ('doublet' if x == 2 else 'singlet')
)
classifications = classifications[['cell_barcode', 'classification']]

#save classifications
os.makedirs(f'results/hashsolo/{dataset_id}', exist_ok=True)
classifications.to_csv(f'results/hashsolo/{dataset_id}/classifications.csv', index=False)

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
summary.to_csv(f'results/hashsolo/{dataset_id}/summary.csv', index=False)

print(summary)
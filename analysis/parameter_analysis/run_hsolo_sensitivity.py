# script to run parameter scan of HashSolo on Mcginnis_MS

##########################################
#      TUNABLE PARAMS:
#        *priors - prior probabilities [negative, singlet, doublet] (default [0.01, 0.8, 0.19]) 
#        *number_of_noise_barcodes - number of barcodes used to create noise distribution. (default len(cell_hashing_columns) - 2)
##########################################

# TO RUN:
#        conda activate demux-py
#        python3 analysis/parameter_analysis/run_hsolo_sensitivity.py

import sys
import anndata as ad
import pandas as pd
import scanpy.external as sce
import itertools
import os

#data loading
input_file = 'data/mcginnis_ms/hto/GSM4904942_8donor_PBMC_AH_MULTI_matrix.csv'
data = pd.read_csv(input_file, index_col=0)
data = data.drop(columns=[col for col in data.columns if 'nUMI' in col])
data = data.iloc[:, :8]

barcode_map = pd.read_csv('config/barcode_maps/mcginnis_ms.csv')
barcode_lookup = dict(zip(barcode_map['barcode'], barcode_map['index']))

adata = ad.AnnData(X=data.values,
                   obs=pd.DataFrame(index=data.index),
                   var=pd.DataFrame(index=data.columns))

for col in data.columns:
    adata.obs[col] = data[col].values

os.makedirs('analysis/parameter_analysis/hashsolo_runs', exist_ok=True)


#parameter scan
prior_neg_vals = [0.01, 0.05, 0.1, 0.15, 0.2]
prior_singlet_vals = [0.7, 0.75, 0.8, 0.85, 0.9]
noise_barcode_vals = [4, 5, 6, 7]  # for 8 barcodes, default is 6

for prior_neg, prior_singlet, n_noise in itertools.product(
    prior_neg_vals, prior_singlet_vals, noise_barcode_vals
):
    prior_doublet = round(1 - prior_neg - prior_singlet, 4)
    if prior_doublet <= 0:
        continue

    adata_copy = adata.copy()
    try:
        sce.pp.hashsolo(
            adata_copy,
            cell_hashing_columns=list(data.columns),
            priors=[prior_neg, prior_singlet, prior_doublet],
            number_of_noise_barcodes=n_noise
        )

        calls = adata_copy.obs['Classification']
        pred = calls.map(lambda x: 0 if x == 'Negative' else (
            1000 if x == 'Doublet' else barcode_lookup.get(x, None)))

        result = pd.DataFrame({
            'cell_barcode': adata_copy.obs.index,
            'assignment': pred,
            'prior_neg': prior_neg,
            'prior_singlet': prior_singlet,
            'prior_doublet': prior_doublet,
            'n_noise': n_noise
        })

        fname = f'analysis/parameter_analysis/hashsolo_runs/hashsolo_{prior_neg}_{prior_singlet}_{n_noise}.csv'
        result.to_csv(fname, index=False)
        print(f'Done: {prior_neg}, {prior_singlet}, {n_noise}')

    except Exception as e:
        print(f'Failed: {prior_neg}, {prior_singlet}, {n_noise} - {e}')
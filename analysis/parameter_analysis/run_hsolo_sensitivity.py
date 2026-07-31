# script to run parameter scan of HashSolo on Mcginnis_MS and Mcginnis_AB

##########################################
#      TUNABLE PARAMS:
#        *priors - prior probabilities [negative, singlet, doublet] (default [0.01, 0.8, 0.19]) 
#        *number_of_noise_barcodes - number of barcodes used to create noise distribution. (default len(cell_hashing_columns) - 2)
##########################################

# TO RUN:
#        conda activate demux-py
#        python3 analysis/parameter_analysis/run_hsolo_sensitivity.py

import anndata as ad
import pandas as pd
import scanpy.external as sce
import itertools
import os

def run_hashsolo_scan(data, barcode_lookup, output_dir,
                      prior_neg_vals, prior_singlet_vals, noise_barcode_vals):
    
    adata = ad.AnnData(X=data.values,
                       obs=pd.DataFrame(index=data.index),
                       var=pd.DataFrame(index=data.columns))
    for col in data.columns:
        adata.obs[col] = data[col].values

    os.makedirs(output_dir, exist_ok=True)

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

            fname = f'{output_dir}/hashsolo_{prior_neg}_{prior_singlet}_{n_noise}.csv'
            result.to_csv(fname, index=False)
            print(f'Done: {prior_neg}, {prior_singlet}, {n_noise}')

        except Exception as e:
            print(f'Failed: {prior_neg}, {prior_singlet}, {n_noise} - {e}')

# parameter grid
prior_neg_vals = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4]
prior_singlet_vals = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95]
noise_barcode_vals = [4, 5, 6, 7]

# ---- mcginnis_ms ----
print('Running mcginnis_ms...')
data_ms = pd.read_csv('data/mcginnis_ms/hto/GSM4904942_8donor_PBMC_AH_MULTI_matrix.csv', index_col=0)
data_ms = data_ms.drop(columns=[col for col in data_ms.columns if 'nUMI' in col])
data_ms = data_ms.iloc[:, :8]

barcode_map_ms = pd.read_csv('config/barcode_maps/mcginnis_ms.csv')
barcode_lookup_ms = dict(zip(barcode_map_ms['barcode'], barcode_map_ms['index']))

run_hashsolo_scan(data_ms, barcode_lookup_ms,
                  'analysis/parameter_analysis/hashsolo_runs/mcginnis_ms',
                  prior_neg_vals, prior_singlet_vals, noise_barcode_vals)

# ---- mcginnis_ab ----
print('Running mcginnis_ab...')
data_ab = pd.read_csv('data/mcginnis_hto/hto/GSM4904939_8donor_PBMC_AH_SCMK_matrix.csv', index_col=0)
data_ab = data_ab.drop(columns=[col for col in data_ab.columns if 'nUMI' in col])
data_ab = data_ab.iloc[:, :8]

barcode_map_ab = pd.read_csv('config/barcode_maps/mcginnis_ab.csv')
barcode_lookup_ab = dict(zip(barcode_map_ab['barcode'], barcode_map_ab['index']))

run_hashsolo_scan(data_ab, barcode_lookup_ab,
                  'analysis/parameter_analysis/hashsolo_runs/mcginnis_ab',
                  prior_neg_vals, prior_singlet_vals, noise_barcode_vals)
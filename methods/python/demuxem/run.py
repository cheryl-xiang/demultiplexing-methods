#script to run DemuxEM

#to run in terminal: 
#    (1) conda activate demux-py
#    (2) python3 methods/python/demuxem/run.py dataset data/dataset/hto/file_name.csv data/dataset/rna/raw_gene_bc_matrices_h5.h5 [switch_transpose]
#    switch_transpose: TRUE to switch default transposing behavior (where barcodes are cols)

import sys
import os
import pandas as pd
import pegasusio as io
from demuxEM import estimate_background_probs, demultiplex

if __name__ == '__main__':
    #read command line arguments
    dataset_id = sys.argv[1]
    input_hto_file = sys.argv[2]
    input_rna_file = sys.argv[3]

    switch_transpose = sys.argv[4].lower() == 'true' if len(sys.argv) >= 5 else False

    # set up output directory and output name prefix
    output_dir = f'results/demuxem/{dataset_id}'
    os.makedirs(output_dir, exist_ok=True)
    output_name = f'{output_dir}/{dataset_id}'

    #data loading
    rna = io.read_input(input_rna_file)
    hto = io.read_input(input_hto_file, transpose=True, modality='hashing')

    if switch_transpose:
        #rna = rna.T
        hto = io.read_input(input_hto_file, transpose=False, modality='hashing')

    print('RNA shape:', rna.shape)
    print('HTO shape:', hto.shape)
    print('Intersection:', len(set(rna.obs_names) & set(hto.obs_names)))

    #run demuxEM
    estimate_background_probs(hto, random_state=0)
    demultiplex(rna, hto, n_threads=1)

    #parse results
    classifications = pd.DataFrame({
        'cell_barcode': rna.obs.index,
        'classification': rna.obs['demux_type'].values
    })

    #standardize labels
    classifications['classification'] = classifications['classification'].map({
        'singlet': 'singlet',
        'doublet': 'doublet',
        'unknown': 'negative'
    })

    #save classifications
    classifications.to_csv(f'{output_dir}/classifications.csv', index=False)

    #save summary counts
    summary = classifications.groupby('classification').size().reset_index(name='n')
    summary['dataset'] = dataset_id
    summary['method'] = 'demuxem'

    total = pd.DataFrame([{
        'classification': 'total',
        'n': summary['n'].sum(),
        'dataset': dataset_id,
        'method': 'demuxem'
    }])

    summary = pd.concat([summary, total], ignore_index=True)
    summary.to_csv(f'{output_dir}/summary.csv', index=False)

    print(summary)
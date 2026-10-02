# convert HTO CSVs to MEX format for cellbouncer

# TO RUN:
#       conda activate demux-py
#       python3 methods/demuxtags/csv_to_mex.py

import pandas as pd
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
import os
import re

datasets = {
    'mcginnis_ms': ('data/mcginnis_ms/hto/GSM4904942_8donor_PBMC_AH_MULTI_matrix.csv', False),
    'mcginnis_ab': ('data/mcginnis_hto/hto/GSM4904939_8donor_PBMC_AH_SCMK_matrix.csv', False),
    'bar11': ('data/ah1/hto/barcodes_matrix_filtered_ah1.csv', False),
    'gaublomme': ('data/gaublomme/hto/barcodes_filtered_gaublomme.csv', True),
    'howitt_batch1_c1': ('data/howitt_b1c1/batch1_c1_hto_counts_howitt.csv', True),
    'howitt_batch1_c2': ('data/howitt_b1c2/batch1_c2_hto_counts_howitt.csv', True),
    'howitt_batch2_c1': ('data/howitt_b2c1/batch2_c1_hto_counts_howitt.csv', True),
    'howitt_batch2_c2': ('data/howitt_b2c2/batch2_c2_hto_counts_howitt.csv', True),
    'howitt_batch3_c1': ('data/howitt_b3c1/batch3_c1_hto_counts_howitt.csv', True),
    'howitt_batch3_c2': ('data/howitt_b3c2/batch3_c2_hto_counts_howitt.csv', True),
    'howitt_capture1': ('data/howitt_cell_cap1/lmo_counts_capture1_howitt_cell_line.csv', True),
    'howitt_capture2': ('data/howitt_cell_cap2/lmo_counts_capture2_howitt_cell_line.csv', True),
    'howitt_capture3': ('data/howitt_cell_cap3/lmo_counts_capture3_howitt_cell_line.csv', True),
    'cook_mix1' : ('data/cook_mix1/GSE147405_TimeCourse_Mix1_barcode_counts_cook.csv', False),
    'cook_mix2' : ('data/cook_mix2/GSE147405_TimeCourse_Mix2_barcode_counts_cook.csv', False),
    'cook_mix3a' : ('data/cook_mix3a/GSE147405_TimeCourse_Mix3a_barcode_counts_cook.csv',False),
    'cook_mix3b' : ('data/cook_mix3b/GSE147405_TimeCourse_Mix3b_barcode_counts_cook.csv',False),
    'cook_mix4a' : ('data/cook_mix4a/GSE147405_TimeCourse_Mix4a_barcode_counts_cook.csv',False),
    'cook_mix4b' : ('data/cook_mix4b/GSE147405_TimeCourse_Mix4b_barcode_counts_cook.csv',False),
}

for dataset_id, (input_file, switch_transpose) in datasets.items():
    print(f'Converting {dataset_id}...')

    data = pd.read_csv(input_file, index_col=0)
    data = data.drop(columns=[col for col in data.columns if 'nUMI' in col])
    data = data.drop(columns=[col for col in data.columns if 'TSNE' in col])

    if switch_transpose:
        data = data.T

    cell_barcodes = data.index.tolist()
    cell_barcodes = [re.sub(r'^\d+_', '', re.sub(r'-\d+$', '', bc)) for bc in cell_barcodes]

    hto_names = data.columns.tolist()

    # MEX format: rows=features(HTOs), cols=cells
    mat = csr_matrix(data.values.T)

    output_dir = f'data/{dataset_id}/hto/mex'
    os.makedirs(output_dir, exist_ok=True)

    with open(f'{output_dir}/barcodes.tsv', 'w') as f:
        for bc in cell_barcodes:
            f.write(bc + '\n')

    with open(f'{output_dir}/features.tsv', 'w') as f:
        for hto in hto_names:
            hto_clean = hto.replace(' ', '_')
            f.write(f'{hto_clean}\t{hto_clean}\tMultiplexing Capture\n')

    mmwrite(f'{output_dir}/matrix.mtx', mat)

    print(f'  {len(cell_barcodes)} cells, {len(hto_names)} HTOs -> {output_dir}')
    print()
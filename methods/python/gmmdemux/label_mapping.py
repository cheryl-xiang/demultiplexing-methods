#script to map gmmdemux assignments and parse runtimes

#TO RUN:
#       conda activate demux-py
#       python3 methods/python/gmmdemux/label_mapping.py

import pandas as pd
import os
import re

datasets = {
    'mcginnis_ms': 8,
    'mcginnis_ab': 8,
    'bar11': 11,
    'gaublomme': 8,
    'howitt_batch1_c1': 8,
    'howitt_batch1_c2': 8,
    'howitt_batch2_c1': 8,
    'howitt_batch2_c2': 8,
    'howitt_batch3_c1': 8,
    'howitt_batch3_c2': 8,
    'howitt_capture1': 3,
    'howitt_capture2': 3,
    'howitt_capture3': 3,
}

for dataset_id, n_htos in datasets.items():
    output_dir = f'results/gmmdemux/{dataset_id}'
    csv_file = os.path.join(output_dir, 'GMM_simplified.csv')

    if not os.path.exists(csv_file):
        print(f'Skipping {dataset_id} - no results found')
        continue

    classifications = pd.read_csv(csv_file, index_col=0)
    classifications.index.name = 'cell_barcode'
    classifications = classifications.reset_index()
    classifications.columns = ['cell_barcode', 'label', 'probability']

    # strip barcode prefixes/suffixes
    classifications['cell_barcode'] = classifications['cell_barcode'].str.replace(r'\.\d+$', '', regex=True)
    classifications['cell_barcode'] = classifications['cell_barcode'].str.replace(r'-\d+$', '', regex=True)
    classifications['cell_barcode'] = classifications['cell_barcode'].str.replace(r'^\d+_', '', regex=True)

    def map_label_string(label, n=n_htos):
        if label == 0 or label == n + 2:
            return 'negative'
        elif label == n + 1:
            return 'doublet'
        elif 1 <= label <= n:
            return 'singlet'
        else:
            return 'negative'

    def map_label_numeric(label, n=n_htos):
        if label == 0 or label == n + 2:
            return 0
        elif label == n + 1:
            return 1000
        elif 1 <= label <= n:
            return label
        else:
            return 0

    classifications['classification'] = classifications['label'].apply(map_label_string)

    classifications[['cell_barcode', 'classification']].to_csv(
        f'{output_dir}/classifications.csv', index=False)

    assignments = classifications[['cell_barcode']].copy()
    assignments['assignment'] = classifications['label'].apply(map_label_numeric)
    assignments.to_csv(f'{output_dir}/assignments.csv', index=False)

    summary = classifications.groupby('classification').size().reset_index(name='n')
    summary['dataset'] = dataset_id
    summary['method'] = 'gmmdemux'

    total = pd.DataFrame([{
        'classification': 'total',
        'n': summary['n'].sum(),
        'dataset': dataset_id,
        'method': 'gmmdemux'
    }])

    summary = pd.concat([summary, total], ignore_index=True)
    summary.to_csv(f'{output_dir}/summary.csv', index=False)

    # parse runtime from time.txt
    time_file = os.path.join(output_dir, 'time.txt')
    if os.path.exists(time_file):
        with open(time_file, 'r') as f:
            content = f.read()

        match = re.search(r'real\s+(\d+)m([\d.]+)s', content)
        if match:
            minutes = int(match.group(1))
            seconds = float(match.group(2))
            runtime = minutes * 60 + seconds
        else:
            match = re.search(r'([\d.]+)\s+total', content)
            if match:
                runtime = float(match.group(1))
            else:
                print(f'{dataset_id}: could not parse time format')
                runtime = None

        if runtime is not None:
            os.makedirs(output_dir, exist_ok=True)
            runtime_df = pd.DataFrame([{
                'dataset': dataset_id,
                'method': 'gmmdemux',
                'runtime_seconds': runtime
            }])
            runtime_df.to_csv(f'{output_dir}/runtime.csv', index=False)
            print(f'{dataset_id} runtime: {runtime:.2f} seconds')
    else:
        print(f'{dataset_id}: no time file found')

    # clean up temp file
    temp_file = f'{output_dir}/temp_hto.csv'
    if os.path.exists(temp_file):
        os.remove(temp_file)

    print(f'{dataset_id}:')
    print(summary)
    print()
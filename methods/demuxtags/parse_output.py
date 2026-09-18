# script to parse demux_tags output

# TO RUN:
#       conda activate demux-py
#       python3 methods/demuxtags/parse_output.py

import pandas as pd
import os
import re

barcode_maps = {
    'mcginnis_ms': {'Bar9': 1, 'Bar10': 2, 'Bar11': 3, 'Bar12': 4, 'Bar13': 5, 'Bar14': 6, 'Bar15': 7, 'Bar16': 8},
    'mcginnis_ab': {'Bar1': 1, 'Bar2': 2, 'Bar3': 3, 'Bar4': 4, 'Bar5': 5, 'Bar6': 6, 'Bar7': 7, 'Bar8': 8},
    'bar11': {'Bar1': 1, 'Bar2': 2, 'Bar3': 3, 'Bar4': 4, 'Bar5': 5, 'Bar6': 6, 'Bar7': 7, 'Bar8': 8, 'Bar9': 9, 'Bar10': 10, 'Bar11': 11},
    'gaublomme': {'S1HuF': 1, 'S2HuM': 2, 'S3HuF': 3, 'S4HuM': 4, 'S5HuF': 5, 'S6HuM': 6, 'S7HuF': 7, 'S8HuM': 8},
    'howitt_batch1_c1': {'BAL_01': 1, 'BAL_02': 2, 'BAL_03': 3, 'BAL_04': 4, 'BAL_05': 5, 'BAL_06': 6, 'BAL_07': 7, 'BAL_08': 8},
    'howitt_batch1_c2': {'BAL_01': 1, 'BAL_02': 2, 'BAL_03': 3, 'BAL_04': 4, 'BAL_05': 5, 'BAL_06': 6, 'BAL_07': 7, 'BAL_08': 8},
    'howitt_batch2_c1': {'BAL_09': 1, 'BAL_10': 2, 'BAL_11': 3, 'BAL_12': 4, 'BAL_13': 5, 'BAL_14': 6, 'BAL_15': 7, 'BAL_16': 8},
    'howitt_batch2_c2': {'BAL_09': 1, 'BAL_10': 2, 'BAL_11': 3, 'BAL_12': 4, 'BAL_13': 5, 'BAL_14': 6, 'BAL_15': 7, 'BAL_16': 8},
    'howitt_batch3_c1': {'BAL_17': 1, 'BAL_18': 2, 'BAL_19': 3, 'BAL_20': 4, 'BAL_21': 5, 'BAL_22': 6, 'BAL_23': 7, 'BAL_24': 8},
    'howitt_batch3_c2': {'BAL_17': 1, 'BAL_18': 2, 'BAL_19': 3, 'BAL_20': 4, 'BAL_21': 5, 'BAL_22': 6, 'BAL_23': 7, 'BAL_24': 8},
    'howitt_capture1': {'CL_01': 1, 'CL_02': 2, 'CL_03': 3},
    'howitt_capture2': {'CL_01': 1, 'CL_02': 2, 'CL_03': 3},
    'howitt_capture3': {'CL_01': 1, 'CL_02': 2, 'CL_03': 3},
}

datasets = list(barcode_maps.keys())

for dataset_id in datasets:
    output_dir = f'results/demuxtags/{dataset_id}'
    assignments_file = f'{output_dir}/{dataset_id}.assignments'
    time_file = f'{output_dir}/time.txt'

    if not os.path.exists(assignments_file):
        print(f'Skipping {dataset_id} - no assignments file found')
        continue

    # load assignments
    df = pd.read_csv(assignments_file, sep='\t', header=None,
                     names=['cell_barcode', 'identity', 'type', 'llr'])

    barcode_lookup = barcode_maps[dataset_id]

    # map to numeric assignments
    def map_assignment(row):
        if row['type'] == 'S':
            return barcode_lookup.get(row['identity'], None)
        elif row['type'] in ['D', 'M']:
            return 1000
        else:
            return 0

    df['assignment'] = df.apply(map_assignment, axis=1)

    # load all barcodes from mex to add missing cells as negative
    barcodes_file = f'data/{dataset_id}/hto/mex/barcodes.tsv'
    all_barcodes = pd.read_csv(barcodes_file, header=None)[0].tolist()
    assigned_barcodes = set(df['cell_barcode'].tolist())
    missing_barcodes = [b for b in all_barcodes if b not in assigned_barcodes]

    if len(missing_barcodes) > 0:
        missing_df = pd.DataFrame({
            'cell_barcode': missing_barcodes,
            'assignment': 0
        })
        assignments_out = pd.concat([df[['cell_barcode', 'assignment']], missing_df], ignore_index=True)
    else:
        assignments_out = df[['cell_barcode', 'assignment']]

    os.makedirs(output_dir, exist_ok=True)
    assignments_out.to_csv(f'{output_dir}/assignments.csv', index=False)

    # save summary
    def map_classification(assignment):
        if assignment == 0:
            return 'negative'
        elif assignment == 1000:
            return 'doublet'
        else:
            return 'singlet'

    assignments_out['classification'] = assignments_out['assignment'].apply(map_classification)
    summary = assignments_out.groupby('classification').size().reset_index(name='n')
    summary['dataset'] = dataset_id
    summary['method'] = 'demuxtags'

    total = pd.DataFrame([{
        'classification': 'total',
        'n': summary['n'].sum(),
        'dataset': dataset_id,
        'method': 'demuxtags'
    }])

    summary = pd.concat([summary, total], ignore_index=True)
    summary.to_csv(f'{output_dir}/summary.csv', index=False)

    # parse runtime
    if os.path.exists(time_file):
        with open(time_file, 'r') as f:
            content = f.read()
        match = re.search(r'real\s+(\d+)m([\d.]+)s', content)
        if match:
            runtime = int(match.group(1)) * 60 + float(match.group(2))
        else:
            match = re.search(r'([\d.]+)\s+total', content)
            runtime = float(match.group(1)) if match else None

        if runtime is not None:
            runtime_df = pd.DataFrame([{
                'dataset': dataset_id,
                'method': 'demuxtags',
                'runtime_seconds': runtime
            }])
            runtime_df.to_csv(f'{output_dir}/runtime.csv', index=False)
            print(f'{dataset_id} runtime: {runtime:.2f} seconds')

    print(f'{dataset_id}:')
    print(summary)
    print()
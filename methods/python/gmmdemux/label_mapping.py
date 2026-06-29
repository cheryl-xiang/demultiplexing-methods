import pandas as pd
import os

datasets = {
    'mcginnis_ms': 8,
    'mcginnis_hto': 8,
    'winkler_pdx1a': 11,
    'winkler_pdx1b': 11,
    'winkler_pdx1c': 11,
    'winkler_pdx1d': 11,
    'stoeckius': 8,
    'ah1': 11,
    'gaublomme': 8,
    'howitt_b1c1': 8,
    'howitt_b1c2': 8,
    'howitt_b2c1': 8,
    'howitt_b2c2': 8,
    'howitt_b3c1': 8,
    'howitt_b3c2': 8,
    'howitt_cell_cap1': 3,
    'howitt_cell_cap2': 3,
    'howitt_cell_cap3': 3,
}

for dataset_id, n_htos in datasets.items():
    output_dir = f'results/gmmdemux/{dataset_id}'
    csv_file = os.path.join(output_dir, 'GMM_simplified.csv')

    if not os.path.exists(csv_file):
        print(f'Skipping {dataset_id} - no results found')
        continue

    if n_htos is None:
        print(f'Skipping {dataset_id} - n_htos not set')
        continue

    classifications = pd.read_csv(csv_file, index_col=0)
    classifications.index.name = 'cell_barcode'
    classifications = classifications.reset_index()
    classifications.columns = ['cell_barcode', 'label', 'probability']

    def map_label(label, n=n_htos):
        if label == 0 or label == n + 2:
            return 'negative'
        elif label == n + 1:
            return 'doublet'
        elif 1 <= label <= n:
            return 'singlet'
        else:
            return 'negative'

    classifications['classification'] = classifications['label'].apply(map_label)
    classifications.to_csv(f'{output_dir}/classifications.csv', index=False)

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

    # clean up temp file
    temp_file = f'{output_dir}/temp_hto.csv'
    if os.path.exists(temp_file):
        os.remove(temp_file)

    print(f'{dataset_id}:')
    print(summary)
    print()
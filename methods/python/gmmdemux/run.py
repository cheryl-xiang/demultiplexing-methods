#script to run GMMDemux

#to run in terminal: 
#    (1) conda activate demux-py
#    (2) python3 methods/python/gmmdemux/run.py dataset data/dataset/hto/file_name.csv

import sys
import os
import subprocess
import pandas as pd
import shutil
import glob

#read command line arguments
dataset_id = sys.argv[1]
input_file = sys.argv[2]

#data loading
data = pd.read_csv(input_file, index_col=0)
data = data.drop(columns=[col for col in data.columns if 'nUMI' in col])
data = data.loc[:, data.sum() > 0]

#get hashtag column names
hto_names = ','.join(data.columns.tolist())

# set up output directory
output_dir = f"results/gmmdemux/{dataset_id}"
os.makedirs(output_dir, exist_ok=True)

#run GMMDemux (from cmd line)
cmd = [
    'GMM-demux',
    '-c', input_file,        
    hto_names,
    '-s', output_dir,
    '-t', '0.8',
    '-rs', '42'
]

subprocess.run(cmd, check=True)

#parse simplified report
csv_file = os.path.join(output_dir, 'GMM_simplified.csv')

classifications = pd.read_csv(csv_file, index_col=0)
classifications.index.name = 'cell_barcode'
classifications = classifications.reset_index()
classifications.columns = ['cell_barcode', 'label', 'probability']

#map simplified labels
def map_label(label):
    if label == 9:
        return 'doublet'
    elif label == 0 or label == 10:
        return 'negative'
    else:
        return 'singlet'

classifications['classification'] = classifications['label'].apply(map_label)

#save classifications
classifications.to_csv(f'{output_dir}/classifications.csv', index=False)

#save summary counts
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
summary.to_csv(f"{output_dir}/summary.csv", index=False)

#delete extra files produced
if os.path.exists('SSD_mtx'):
    shutil.rmtree('SSD_mtx')

for folder in glob.glob('GMM_Demux_*'):
    shutil.rmtree(folder)

print(summary)
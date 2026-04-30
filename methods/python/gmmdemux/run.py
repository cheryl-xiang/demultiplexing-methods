#script to run GMMDemux

#to run in terminal: 
#    (1) conda activate demux-py
#    (2) python3 methods/python/gmmdemux/run.py dataset_# data/file_name.csv

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
data = data.drop(columns=[col for col in data.columns if 'nUMI' in col])  #also check for other col names in other data

#get hashtag column names
hto_names = ','.join(data.columns.tolist())

# set up output directory
output_dir = f"results/gmm_demux/{dataset_id}"
full_report_dir = f"{output_dir}/full_report"
os.makedirs(output_dir, exist_ok=True)

#run GMMDemux (from cmd line)
cmd = [
    'GMM-demux',
    '-c', input_file,        
    hto_names,
    '-f', full_report_dir,  #full report to parse 
    '-s', output_dir        #simplified report
]

subprocess.run(cmd, check=True)

#parse report 
config_file = os.path.join(full_report_dir, 'GMM_full.config')
csv_file = os.path.join(full_report_dir, 'GMM_full.csv')

#map nums back to labels
label_map = {}
with open(config_file, 'r') as f:
    for line in f:
        parts = line.strip().split(', ')
        if len(parts) == 2:
            label_map[parts[0]] = parts[1] 

# read classifications
classifications = pd.read_csv(csv_file, index_col=0)
classifications.index.name = 'cell_barcode'
classifications = classifications.reset_index()
classifications.columns = ['cell_barcode', 'label', 'probability']


#match classfication names
def map_label(label):
    name = label_map.get(str(int(label)), 'negative')
    if name == 'negative':
        return 'negative'
    elif '-' in name:
        return 'multiplet'
    else:
        return 'singlet'
    
classifications['classification'] = classifications['label'].apply(map_label)

#save classifications
classifications.to_csv(f'{output_dir}/classifications.csv', index=False)

#save summary counts
summary = classifications.groupby('classification').size().reset_index(name='n')
summary['dataset'] = dataset_id
summary['method'] = 'gmm_demux'

total = pd.DataFrame([{
    'classification': 'total',
    'n': summary['n'].sum(),
    'dataset': dataset_id,
    'method': 'gmm_demux'
}])

summary = pd.concat([summary, total], ignore_index=True)
summary.to_csv(f"{output_dir}/summary.csv", index=False)

#delete extra files produced
if os.path.exists('SSD_mtx'):
    shutil.rmtree('SSD_mtx')

for folder in glob.glob('GMM_Demux_*'):
    shutil.rmtree(folder)

if os.path.exists(full_report_dir):
    shutil.rmtree(full_report_dir)

if os.path.exists(f'{output_dir}/GMM_simplified.csv'):
    os.remove(f'{output_dir}/GMM_simplified.csv')

print(summary)

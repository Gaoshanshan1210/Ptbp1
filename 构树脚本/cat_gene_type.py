import pandas as pd
import sys
import os

# cat files base on gene type 
dir_i = sys.argv[1] # Directory of input files, main_folder - subfolders - files
dir_o = sys.argv[2] # output directory

gene_types = {}
for sub_folder in os.listdir(dir_i):
    if sub_folder != '.DS_Store':
        print(f'[*] Folder: {sub_folder}')
        curr_folder = os.path.join(dir_i, sub_folder)

        for file in os.listdir(curr_folder):
            if file != '.DS_Store':
                print(f'[*] \tL File: {file}')
                gene_types.setdefault(file.split('.')[1], '')
                curr_file = os.path.join(curr_folder, file)

                try:
                    if type(gene_types[file.split('.')[1]]) != pd.DataFrame:
                        gene_types[file.split('.')[1]] = pd.read_csv(curr_file, sep='\t')

                    else:
                        gene_types[file.split('.')[1]] = pd.concat([gene_types[file.split('.')[1]], pd.read_csv(curr_file, sep='\t')], ignore_index=True)

                except pd.errors.EmptyDataError:
                    pass

for k, v in gene_types.items():
    print(f'[*] Generating output file for {k}...')
    try:
        v.to_csv(os.path.join(dir_o, f'{k}.csv'), sep='\t', index=False)
    except AttributeError:
        pass

print('[*] Done.')

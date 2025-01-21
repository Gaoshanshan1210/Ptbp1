import pandas as pd
import argparse
import sys
import os

class Tree_data_prepare:
    def __init__(self, fi, folder_o, sep, name, ref, min_repeat, max_repeat, recursive):
        self.fi = fi
        self.folder_o = folder_o
        self.sep = sep
        self.name = name
        self.ref = ref
        self.min_repeat = min_repeat
        self.max_repeat = max_repeat
        self.recursive = recursive

    # Update data by reference file
    def convert2Ref(self):
        print('[*] Reading from reference file...')
        ref_df = pd.read_csv(self.ref, sep=self.sep)
        
        # Dictionary for old - new
        update_catelog = {}
        for i in ref_df.index:
            try:
                if ref_df.loc[i, 'IDs'] == self.name:
                    key = ref_df.loc[i, 'key']
                    update_catelog.setdefault(key, '')
                    update_catelog[key] = ref_df.loc[i, 'value']
            
            except KeyError as e:
                sys.exit(f'[*] Error: {e}. Please change the column name of reference file as required, check --help to learn more.')
        print('[*] Done.')

        # Update data
        print('[*] Updating old data to new data...')
        df = pd.read_csv(self.fi, sep=self.sep)
        statistic = {}
        try:
            for i in df.index:
                key = df.loc[i, 'key']
                df.loc[i, 'old'] = update_catelog[key]
                
                statistic.setdefault(key, [])
                statistic[key].append(i)
        
        except KeyError as e:
                sys.exit(f'[*] Error: {e}. Please change the column name of input file as required, check --help to learn more.')
        print('[*] Done.')

        # Filter
        if self.min_repeat != 0 or self.max_repeat != 0:
            print('[*] Performing filting...')
            drop_index = []
            for v in statistic.values():
                if len(v) < self.min_repeat:
                    drop_index += v
                
                if self.max_repeat != 0 and len(v) > self.max_repeat:
                    drop_index += v
            
            df.drop(index=drop_index, inplace=True)
            print('[*] Done.')

        # Save updated file as ID.csv
        print('[*] Saving updated file...')
        if self.recursive == 'T':
            fo_path = os.path.join(self.folder_o, f"{os.path.basename(self.fi).replace('.csv', '').replace('.tsv', '')}.{self.name}.csv")
        else:
            fo_path = os.path.join(self.folder_o, f'{self.name}.csv')
        df.to_csv(fo_path, sep=self.sep, index=False)
        print('[*] Done.')

        return fo_path
    
    # Generate 0,1 matrix
    def to2Dmatrix(self, fi):
        print('[*] Reading file...')
        df = pd.read_csv(fi, sep=self.sep)
        print('[*] Done.')
        
        # Column name: CB, (#CHROM)_(Start)_(REF)>(ALT_expected)...
        # Row name: (CB)_(Cell_type_observed)...
        raw_df = {'CB': []}
        for i in range(df.shape[0]):
            print(f'\r[*] Converting file to 2D matrix: {round(round((i+1)/df.shape[0], 3)*100, 2)}%', end='')

            curr_CB = f'{df["old"][i].replace(" ", "_")}_{df["key"][i]}'
            column_name = f'{df["#CHROM"][i]}_{df["Start"][i]}_{df["REF"][i]}>{df["ALT_expected"][i]}'

            # CB must not repeated
            if curr_CB not in raw_df['CB']:
                raw_df['CB'].append(curr_CB)

                # For current column
                raw_df.setdefault(column_name, [])
                
                # keep the length same as column CB
                completement = len(raw_df['CB']) - len(raw_df[column_name])
                if completement > 0:
                    for j in range(completement - 1):
                        raw_df[column_name].append(0)

                raw_df[column_name].append(1)

                # For other columns, fill 0
                for column in raw_df.keys():
                    if column != 'CB' and column != column_name:
                        # keep the length same as column CB
                        completement = len(raw_df['CB']) - len(raw_df[column])
                        if completement > 0:
                            for j in range(completement):
                                raw_df[column].append(0)
            
            else:
                index = raw_df['CB'].index(curr_CB)
                
                raw_df.setdefault(column_name, [])
                completement = len(raw_df['CB']) - len(raw_df[column_name])
                if completement > 0:
                    for j in range(completement):
                        raw_df[column_name].append(0)

                raw_df[column_name][index] = 1
        
        print('\n[*] Done.')

        # length check
        #for k, v in raw_df.items():
        #    print(f'{k}: {len(v)}')

        # Generating output
        print('[*] Generating output file...')
        try:
            # Add 7 non-mutate with all 0
            for k in raw_df.keys():
                if k == 'CB':
                    for i in range(8):
                        #print(f'{column}: {len(raw_df[column])}')
                        raw_df[k].append(f'non-mutated_{i}')
                        #print(f'L {column}: {len(raw_df[column])}')

                else:
                    for i in range(8):
                        #print(f'{column}: {len(raw_df[column])}')
                        raw_df[k].append(0)
                        #print(f'L {column}: {len(raw_df[column])}')

            new_df = pd.DataFrame(raw_df)
        
        except Exception as E:
            sys.exit(f'[*] Error: {E}.')
        
        else:
            if self.recursive == 'T':
                fo_path = os.path.join(self.folder_o, f"matrix_{os.path.basename(self.fi).replace('.csv', '').replace('.tsv', '')}.{self.name}.csv")

            else:
                fo_path = os.path.join(self.folder_o, f'matrix_{self.name}.csv')
            new_df.to_csv(fo_path, sep=sep, index=False)
            print('[*] Done.')
        
            return fo_path
    
    # Generate fasta file base on 0,1 matrix
    def matrix2fasta(self, fi):
        print('[*] Gnerating fasta file...')
        df = pd.read_csv(fi, sep=self.sep)

        if self.recursive == 'T':
            fo_path = os.path.join(self.folder_o, f"{os.path.basename(self.fi).replace('.csv', '').replace('.tsv', '')}.{self.name}.fasta")
        else:
            fo_path = os.path.join(self.folder_o, f'{self.name}.fasta')
        with open(fo_path, 'w') as filo:
            for i in range(df.shape[0]):
                filo.write(f'>{df["CB"][i]}\n')
                
                tmp = ''
                for column in df.columns:
                    if column != 'CB':
                        if int(df[column][i]) > 0:
                            tmp += column.split('_')[-1].split('>')[1]
                        else:
                            tmp += column.split('_')[-1].split('>')[0]
                
                filo.write(f'{tmp}\n')
        
        print('[*] Done.')
    
    def run(self, fi):
        if fi != None:
            self.fi = fi

        if self.ref != 'skip':
            fo_path = self.convert2Ref()
        
        else:
            fo_path = self.fi
        
        loop = True
        while loop:
            if self.ref == 'skip':
                ifContinue = 'y'
            else:
                ifContinue = input('[*] Do you wish to continue(y/n): ')

            if ifContinue.lower() == 'y':
                fo_path = self.to2Dmatrix(fo_path)
                self.matrix2fasta(fo_path)
                loop = False
            
            elif ifContinue.lower() == 'n':
                sys.exit('[*] Done.')
        
    def multi_fi(self):
        if self.recursive == 'T':
            dir = self.fi
            for fi in os.listdir(dir):
                print(f'[*] Current file - {fi}')
                curr_fi = os.path.join(dir, fi)

                self.run(curr_fi)

        else:
            self.run(None)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-I', '--INPUT_FILE', required=True, help='Enter the path of input file. Please change the column name as: target column - key(usually "CB")(same as -REF) | old value column - old(usually "Cell_type_observed")')
    parser.add_argument('-O', '--OUTPUT_DIRECTORY', required=True, help='Enter the directory for outputs.')
    parser.add_argument('-SEP', '--SEPARATION', default='\t', help='Enter the separation for all files, all files must use same separation! Default is TAB.')
    parser.add_argument('-NAME', '--INDIVIDUAL_NAME', required=True, help='Enter the name of individual you are going to update, same as the name in main column in reference file.')
    parser.add_argument('-REF', '--REFERENCE_FILE', default='skip', help='Type "skip" to skip update or enter the path of update reference file. Please change the column name as: main column - IDs | target column - key | new value column - value')
    parser.add_argument('-MIN', '--MINIMUM_DUPLICATE', type=int, default=0, help='Set -MIN to fliter rows for -I.')
    parser.add_argument('-MAX', '--MAXIMUM_DUPLICATE', type=int, default=0, help='Set -MAX to fliter rows for -I.')
    parser.add_argument('-R', '--RECURSIVE', default='F', choices=['T', 'F'], help='To process multiple files in a row, use -R T, -I should be a directory path.')

    args = parser.parse_args()
    fi = args.INPUT_FILE
    folder_o = args.OUTPUT_DIRECTORY
    sep = args.SEPARATION
    name = args.INDIVIDUAL_NAME
    ref = args.REFERENCE_FILE
    min_repeat = args.MINIMUM_DUPLICATE
    max_repeat = args.MAXIMUM_DUPLICATE
    recursive = args.RECURSIVE

    if not os.path.isdir(folder_o):
        sys.exit('[*] -O must be a directory, not a file!')
    
    if recursive == 'T':
        if not os.path.isdir(fi):
            sys.exit('[*] As you set -R T, -I must be a directory, not a file!')

    tdp = Tree_data_prepare(fi, folder_o, sep, name, ref, min_repeat, max_repeat, recursive)
    tdp.multi_fi()

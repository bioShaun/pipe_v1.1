'''
Usage:
    merge_diff_table_by_name.py <diff_table_dir> <merged_diff_table>

'''

import pandas as pd
from os import path
from os import listdir
from docopt import docopt 

arguments = docopt(__doc__, version = 'v1')

table_dir = arguments['<diff_table_dir>']
merged_table = arguments['<merged_diff_table>']

MERGE_NUM = 9
HEADER = ['logFC', 'PValue', 'FDR']
COUNT = 0

for each_dir in listdir(table_dir):
    each_dir_path = path.join(table_dir, each_dir)
    if path.isdir(each_dir_path) and COUNT < 9:
        for each_file in listdir(each_dir_path):
            if 'DE_results' not in each_file:
                continue
            each_file_path = path.join(each_dir_path, each_file)
            each_file_df = pd.read_table(each_file_path, sep = '\t', index_col = 0)
            each_file_out_df = each_file_df.loc[:,HEADER]
            each_file_out_df.loc[:,'Compare'] = each_dir
            if COUNT == 0 :
                global merged_df 
                merged_df = each_file_out_df
                COUNT += 1
                continue
            else:
                merged_df = pd.concat([merged_df, each_file_out_df])
                COUNT += 1

merged_df.to_csv(merged_table, sep = '\t', float_format='%.3e')

'''
extract table from id_file

'''

import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--id', help = 'id file', required = True)
parser.add_argument('--table', help = 'table', required = True)
parser.add_argument('--output', help = 'merged table', required = True)
parser.add_argument('--sep', help = 'sep of two tables ', default = '\t')
parser.add_argument('--header', help = 'header of two tables', default = 'no')
parser.add_argument('--id_col', help = 'id column', default = '1')
parser.add_argument('--choice', help = 'in or out', default = 'in')
args = parser.parse_args()

id_file = args.id
table = args.table
output = args.output
sep_stat = args.sep
header_stat = args.header
col_stat = args.id_col
choice = args.choice

id_dict = {}
with open(id_file,'r') as id_info :
    for n,eachline in enumerate(id_info) :
        eachline = eachline.strip()
        id_dict[eachline] = 1

output_file = open(output,'w')
with open(table,'r') as table_info :
    for n,eachline in enumerate(table_info) :
        if header_stat == 'yes' and n == 0 :
            output_file.write(eachline)
        else :
            eachline_info = eachline.strip().split(sep_stat)
            id_col = int(col_stat) - 1
            if eachline_info[id_col] in id_dict and args.choice == 'in':
                output_file.write(eachline)
            elif eachline_info[id_col] not in id_dict and args.choice == 'out':
                output_file.write(eachline)
            else:
                pass
output_file.close()

"""Usage: reorder_table.py -t FILE -s ORDER_LIST -o OUTFILE

Process FILE and reorder FILE according to ORDER_LIST

Arguments:
  FILE        input file
  ORDER_LIST  selected columns 
  OUTFILE     reordered file

Options:
  -h --help
  -t <file>, --table <file>  input table
  -s <list>, --list <list>   reorder table according to the list
  -o <out>, --output <out>   output file  

"""
from docopt import docopt
import pandas as pd
from os import path

if __name__ == '__main__':
    arguments = docopt(__doc__, version = "v1")
    input_table = arguments['--table']
    input_list = arguments['--list']
    output_table = arguments['--output']

    input_df = pd.read_table(input_table, index_col = 0)
    input_list_inf = [each.strip() for each in open(input_list)]
    reorder_df = input_df[input_list_inf]
    reorder_df.to_csv(output_table, sep = "\t", float_format='%.3f')

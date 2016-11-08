#! /usr/bin/python

'''
data:2016-3-21

For kallisto quant and differential analysis

last update: 2016-3-28
'''

import os
import sys
import argparse
import RNAseq_tools

parser = argparse.ArgumentParser()
parser.add_argument('--transcript', help = 'transcript fasta file.',default = '')
parser.add_argument('--index', help = 'transcript fasta file index.', required =True)
parser.add_argument('--fq_type', help = 'fq file directory.', default = 'PE')
parser.add_argument('--fq_dir', help = 'fq file directory.', required = True)
parser.add_argument('--name_map', help = 'file map sample name to fq file.', default = None)
parser.add_argument('--out_dir', help = 'output directory of quant results.', required = True)
args = parser.parse_args()

my_kallisto = RNAseq_tools.run_kallisto()
my_kallisto.index = args.index
my_kallisto.transcript_fa = args.transcript
my_kallisto.out_dir = args.out_dir

my_fq_files = RNAseq_tools.treat_fq()
my_fq_files.fq_flag = args.fq_type
my_fq_files.clean_dir = args.fq_dir
if args.name_map :
    my_fq_files.mapfile = args.name_map
    my_fq_files.lib_name_to_sample_name()
    my_fq_files.get_sample_name_order()
my_fq_files.store_clean_fq()

if __name__ == '__main__':    
    tr_kallisto_index = args.index
    if not os.path.isfile(tr_kallisto_index) :
        if not os.path.isfile(args.transcript) :
            sys.exit('you need to have either transcript fasta or transcript index')
        tr_kallisto_index = my_kallisto.kallisto_index()

    if my_fq_files.sample_name_order :
        sample_name_order = my_fq_files.sample_name_order
    else :
        sample_name_order = sorted(my_fq_files.fq_dict.keys())
    for each_sample in sample_name_order :
        each_sample_quant = RNAseq_tools.run_kallisto()
        each_sample_quant.index = tr_kallisto_index
        each_sample_quant.out_dir = os.path.join(args.out_dir,each_sample)
        if args.fq_type == 'SE' :
            each_sample_quant.fq_list = [my_fq_files.fq_dict[each_sample]]
            each_sample_quant.fq_length = my_fq_files.get_fq_len(my_fq_files.fq_dict[each_sample])            
        elif args.fq_type == 'PE' :
            each_sample_quant.fq_list = [my_fq_files.fq_dict[each_sample]['1'],my_fq_files.fq_dict[each_sample]['2']]
        else :
            sys.exit('wrong fq type! : %s' % fq_type )
        each_sample_quant.kallisto_quant()

    my_kallisto_result = RNAseq_tools.run_kallisto()
    my_kallisto_result.out_dir = args.out_dir
    my_kallisto_result.sample_list = sample_name_order
    my_kallisto_result.kallisto_merge()










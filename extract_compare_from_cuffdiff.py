#! /usr/bin/python

'''
data:2016-3-14

1. Extract each compare diff analysis results from cuffdiff output

'''

import sys
import os
import argparse

mycwd = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument('--cuffdiff_out',help = 'cuffdiff diff analysis out_dir file.', required = True)
parser.add_argument('--out_dir', help = 'Compare name list.', default = mycwd)
args = parser.parse_args()

diff_file_name = os.path.basename(args.cuffdiff_out)
diff_tag = diff_file_name.split('_')[0]

diff_exp_dict = {}
for n,eachline in enumerate(open(args.cuffdiff_out)) :
    if n != 0 :
        eachline_info = eachline.strip().split('\t')
        compare_name = '%s_vs_%s' % (eachline_info[5],eachline_info[4])
        test_id = eachline_info[0]
        gene_id = eachline_info[1]
        gene_symbol = eachline_info[2]
        gene_location = eachline_info[3]
        control_fpkm = eachline_info[7]
        treat_fpkm = eachline_info[8]
        l2fc = eachline_info[9]
        pvalue = eachline_info[11]
        qvalue = eachline_info[12]
        exp_line_list = [test_id,gene_id,gene_symbol,gene_location,treat_fpkm,control_fpkm,l2fc,pvalue,qvalue]
        if diff_tag == 'gene' :
            exp_line_list.pop(0)
        output_line = '\t'.join(exp_line_list)
        if compare_name not in diff_exp_dict :
            diff_exp_dict[compare_name] = 1
            compare_diff_file = os.path.join(args.out_dir,'%s.cuffdiff.DE_results' % compare_name)
            with open(compare_diff_file,'w') as compare_diff_out :
                header = "Gene_ID\tGene_Symbol\tGene_location\t%s_fpkm\t%s_fpkm\tlog2foldchange\tpvalue\tqvaule\n" % (eachline_info[5],eachline_info[4])
                if diff_tag != 'gene' :
                    header = 'Transcript_ID\t%s' % header
                compare_diff_out.write(header)
        with open(compare_diff_file,'a') as compare_diff_out :
            compare_diff_out.write('%s\n' % output_line)





import sys
import os
import argparse
import re
import python_tools

parser = argparse.ArgumentParser()
parser.add_argument('--maf_fasta', help = 'maf file.', required = True)
#parser.add_argument('--species', help = 'maf species.', required = True)
parser.add_argument('--out_dir', help = 'output directory', required = True)
args = parser.parse_args()

out_dir = os.path.abspath(args.out_dir)

if not os.path.exists(out_dir):
    python_tools.circ_mkdir_unix(out_dir)

id_dict = {}

outfile_list = os.path.join(out_dir,'maf_fasta.list')
outfile_list_info = open(outfile_list,'w')

with open(args.maf_fasta) as maf_fasta_info:
    for eachline in maf_fasta_info:
        eachline = eachline.strip()
        if '>' in eachline:            
            header = re.sub('.TU','|TU',eachline)
            tr_id = header.split('|')[1]
            if tr_id not in id_dict:
                id_dict[tr_id] = 1
                each_tr_out_file = os.path.join(out_dir,'%s.fa' % tr_id)
                outfile_list_info.write('%s\n' % each_tr_out_file)
                each_tr_out_info = open(each_tr_out_file,'w')
            each_tr_out_info.write('%s\n' % header)
        elif eachline == '':
            each_tr_out_info.close()
        else:
            each_tr_out_info.write('%s\n' % eachline)
outfile_list_info.close()

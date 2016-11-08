#coding=utf-8

'''

Treat large RNAseq data

'''

import argparse
import RNAseq_tools
import python_tools
import os
import sys
import time
import json
import re
import random

cwd = os.getcwd()
parser = argparse.ArgumentParser()
parser.add_argument('--datasize_file', help = 'project datasize file', required = True)
parser.add_argument('--data_dir', help = 'original data directory', required = True)
parser.add_argument('--cutoff', help = 'output data directory', default = '10')
parser.add_argument('--target', help = 'target datasize', default = '6')
args = parser.parse_args()

sample_list = []
cutoff = int(args.cutoff)
target_size = int(args.target)
target_reads_low = target_size*4*(10**9)/300
target_reads_up = (target_size+1)*4*(10**9)/300
target_reads_range = xrange(target_reads_low,target_reads_up)

with open(args.datasize_file) as datasize_file_info:
    for n,eachline in enumerate(datasize_file_info):
        if n!= 0:
            eachline_info = eachline.strip().split('\t')
            sample_id = eachline_info[0]
            datasize = float(eachline_info[1])
            if datasize > cutoff:
                sample_list.append(sample_id)

extract_cmd_file = os.path.join(cwd,'extract_cmd.sh')
extract_cmd_list = []
tmp_dir = os.path.join(args.data_dir,'tmp')
if not os.path.exists(tmp_dir):
    extract_cmd_list.append('mkdir -p "%s"' % tmp_dir)

sample_datasize = random.sample(target_reads_range,len(sample_list))
for n,each_sample in enumerate(sample_list):
    extract_cmd_list.append('echo "#### start extract data of sample %s ####"' % each_sample)
    extract_cmd_list.append('date')
    for read_num in (1,2):
        each_fq_file = os.path.join(args.data_dir,'%s_R%s.fastq.gz' % (each_sample,read_num))
        each_fq_bak_file = os.path.join(tmp_dir,'%s_R%s.fastq.gz' % (each_sample,read_num))
        each_extract_read = sample_datasize[n]
        extract_cmd_list.append('mv "%s" "%s"' % (each_fq_file,each_fq_bak_file))        
        extract_cmd_list.append('zcat "%s" | head -%s |gzip > "%s"' % (each_fq_bak_file,each_extract_read,each_fq_file))
    extract_cmd_list.append('date')
    extract_cmd_list.append('echo "#### finish extract data of sample %s ####"\n' % each_sample)
python_tools.write_obj_to_file(extract_cmd_list,extract_cmd_file)



#coding=utf-8

'''

cp RNAseq and Resequencing data to custumer

'''

import argparse
import RNAseq_tools
import python_tools
import os
import sys
import time
import json
import re
import python_tools

cwd = os.getcwd()
parser = argparse.ArgumentParser()
parser.add_argument('--proj_names', help = 'project names, seperate with ","', required = True)
parser.add_argument('--out_dir', help = 'output directory', required = True)
parser.add_argument('--seq_num',help = 'data seq numbers ,seperated with ","', default = "")
parser.add_argument('--wgc_id',help = 'wgc id file ', default = "")
parser.add_argument('--md5_id',help = 'md5 id file ', default = "")
parser.add_argument('--ALL', action='store_true', help = 'cp data in all seq number', default = False)
parser.add_argument('--type', help = 'A rnaseq project or a resequencing project',choices = ['rna','dna'], required = True)
args = parser.parse_args()


data_dir_list_file = '/home/public/project/rawdata/statistics/rawdata_dir_list'
# data_dir_list = python_tools.file_to_list(data_dir_list_file)
wgc_id_list = []
if args.wgc_id:
    wgc_id_list = [each.strip() for each in open(args.wgc_id)] 
if args.md5_id:
    md5_id_list = [each.strip() for each in open(args.md5_id)]

data_dir_dict = {}
with open(data_dir_list_file) as data_dir_list_file_info:
    for eachline in data_dir_list_file_info:
        eachline_info = eachline.strip().split('\t')
        seq_num = eachline_info[0]
        seq_dir = eachline_info[1]
        python_tools.add_dict_value(data_dir_dict,seq_num,seq_dir)
seq_number_list = args.seq_num.split(',')

wgc_to_sample_dict = {}
sample_to_wgc_dict = {}
all_info_file = '/home/lxgui/Documents/paiji/check_up_dowm/出库样品信息.txt'

proj_name_list = args.proj_names.split(',')

with open(all_info_file) as all_info_file_if :
    for n,eachline in enumerate(all_info_file_if):
        if n != 0:
            eachline_info = eachline.strip().split('\t')
            if eachline_info[0] in proj_name_list:
                wgc_id = eachline_info[1]
                sample_id = eachline_info[4]
                wgc_to_sample_dict[wgc_id] = sample_id
                sample_to_wgc_dict[sample_id] = wgc_id

if os.path.isfile('rawdata_number.json') and os.stat('rawdata_number.json').st_size :
    sample_number_dict = python_tools.load_fn_to_obj('rawdata_number.json')
else :
    sample_number_dict = {}

sample_dict = {}
for each_num in data_dir_dict:
    if args.ALL or each_num in seq_number_list:
        each_dir_list = data_dir_dict[each_num]
        for each_dir in each_dir_list:
            rawdata_list = os.listdir(each_dir)
            for each_rawdata_dir in rawdata_list:
                each_fq_dir = os.path.join(each_dir,each_rawdata_dir)
                if os.path.isdir(each_fq_dir):            
                    if re.search(r'(WGC\d{6})',each_rawdata_dir):
                        wgc_id = re.search(r'(WGC\d{6})',each_rawdata_dir).groups()[0]
                        if wgc_id in wgc_to_sample_dict:
                            sample_id = wgc_to_sample_dict[wgc_id]
                            if wgc_id_list and wgc_id not in wgc_id_list:
                                continue
                        else:
                            continue
                    else:
                        sample_id = each_rawdata_dir
                        if sample_id not in sample_to_wgc_dict:
                            continue
                    fq_file_list = os.listdir(each_fq_dir)
                    for each_fq in fq_file_list:
                        each_fq_path = os.path.join(each_fq_dir,each_fq)
                        if sample_id not in sample_dict:
                            sample_dict[sample_id] = RNAseq_tools.rawdata()
                            sample_dict[sample_id].name = sample_id
                        else :
                            pass
                        if each_fq_path.endswith('R1.fastq.gz') and each_fq_path not in sample_dict[sample_id].read1:
                            sample_dict[sample_id].read1.append(each_fq_path)
                        elif each_fq_path.endswith('R2.fastq.gz') and each_fq_path not in sample_dict[sample_id].read2:
                            sample_dict[sample_id].read2.append(each_fq_path)
                        else :
                            pass

python_tools.circ_mkdir_unix(args.out_dir)
time_info = time.localtime()
output_time = '%s-%s-%s-%s:%s:%s' % (time_info.tm_year,time_info.tm_mon,time_info.tm_mday,time_info.tm_hour,time_info.tm_min,time_info.tm_sec)

cp_cmd = os.path.join(cwd,'%s_cp.sh' % output_time)
cp_data_info_file = os.path.join(cwd,'%s_rawdata.info' % output_time)

total_size = []
cp_data_info = open(cp_data_info_file,'w')
for each in sample_dict :
    if each in sample_number_dict :
        sample_dict[each].pre_num = sample_number_dict[each]
    else :
        sample_number_dict[each] = len(sample_dict[each].read1)
    if args.type == 'rna' :
        cmd_line = sample_dict[each].merge_rawdata(args.out_dir)
    else :
        cmd_line = sample_dict[each].cp_rawdata_to_custum(args.out_dir)
    for n,each_fq in enumerate(sample_dict[each].read1) :
        read1_fq = sample_dict[each].read1[n]
        read2_fq = sample_dict[each].read2[n]
        read1_fq_size = os.stat(read1_fq).st_size/float(1024**3)
        read2_fq_size = os.stat(read2_fq).st_size/float(1024**3)
        total_size.extend([read1_fq_size,read2_fq_size])
        read1_fq_size_out = round(read1_fq_size,2)
        read2_fq_size_out = round(read2_fq_size,2)
        cp_data_info.write('%s\t%s\t%sG\t%s\t%sG\n' % (sample_dict[each].name,read1_fq,read1_fq_size_out,read2_fq,read2_fq_size_out))
    python_tools.write_obj_to_file(cmd_line,cp_cmd,True)
cp_data_info.write('total : %sG' % round(sum(total_size),2))
cp_data_info.close()

cp_data_info_json = os.path.join(cwd,'rawdata_number.json')
python_tools.write_obj_to_json(sample_number_dict,cp_data_info_json)

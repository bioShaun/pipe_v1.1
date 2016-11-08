#coding=utf-8

'''

check Resequencing data md5

'''

import argparse
import RNAseq_tools
import python_tools
import os
import sys
import time
import json
import re
import logging

cwd = os.getcwd()
parser = argparse.ArgumentParser()
parser.add_argument('--proj_names', help = 'project names, seperate with ","', required = True)
parser.add_argument('--out_dir', help = 'output directory', required = True)
parser.add_argument('--seq_num', help = 'seq number', default = "")
parser.add_argument('--type', help = 'A rnaseq project or a resequencing project',choices = ['rna','dna'], required = True)
parser.add_argument('--nocheck',  help = 'not check, just cp md5 information from sequencing directory.', action='store_true', default = False)
parser.add_argument('--ALL', action='store_true', help = 'check/get all data md5', default = False)
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG,
        format='%(asctime)s %(filename)s [line:%(lineno)d] %(levelname)s : %(message)s',
        datefmt='%a, %d %b %Y %H:%M:%S',
        filename='{cwd}/check_md5.log'.format(**locals()),
        filemode='w')

data_dir_list_file = '/home/public/project/rawdata/statistics/rawdata_dir_list'
# data_dir_list = python_tools.file_to_list(data_dir_list_file)
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

md5_file = os.path.join(args.out_dir,'fq_md5.txt')
md5_list = []
for each in sample_dict :
    if args.nocheck:
        logging.info('get md5 of %s start' % each)
        log_list = sample_dict[each].get_dna_md5()
        logging.info('get md5 of %s finished' % each)        
    else:
        logging.info('check md5 of %s start' % each)
        if args.type == 'dna' :
            log_list = sample_dict[each].check_dna_md5(args.out_dir)
        else:
            log_list = sample_dict[each].check_rna_md5(args.out_dir)
    if log_list:
        for each_log in log_list:
            if 'ok' in each_log:
                logging.info(each_log)
            else:
                logging.error(each_log)
        logging.info('check md5 of %s finished' % each)
    md5_list.extend(sample_dict[each].md5)
python_tools.write_obj_to_file(md5_list,md5_file)





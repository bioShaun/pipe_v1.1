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

cwd = os.getcwd()
parser = argparse.ArgumentParser()
parser.add_argument('--sample_map', help = 'id map wgc id to sample id', required = True)
parser.add_argument('--out_dir', help = 'output directory', required = True)
parser.add_argument('--rawdata_list', help = 'directory list file', required = True)
args = parser.parse_args()

python_tools.circ_mkdir_unix(args.out_dir)
time_info = time.localtime()
output_time = '%s-%s-%s-%s:%s:%s' % (time_info.tm_year,time_info.tm_mon,time_info.tm_mday,time_info.tm_hour,time_info.tm_min,time_info.tm_sec)

data_dir_list_file = args.rawdata_list
data_dir_list = python_tools.file_to_list(data_dir_list_file)

wgc_to_sample_dict = python_tools.table_to_dict(args.sample_map,1,2,False,'\t')
sample_to_wgc_dict = python_tools.table_to_dict(args.sample_map,2,1,False,'\t')

cp_cmd = os.path.join(cwd,'%s_cp.sh' % output_time)
cp_data_info_file = os.path.join(cwd,'%s_rawdata.info' % output_time)

def get_fq_sample_name(fqname) :
    if '_R1.fastq.gz' in fqname :
        return 'R1',fqname.split('_clean_R1.fastq.gz')[0]
    else :
        return 'R2',fqname.split('_clean_R2.fastq.gz')[0]

first_rawdata_list = os.listdir(args.out_dir)

prepare_sh = open('%s_pre.sh' % output_time,'w')
for each in first_rawdata_list :
    old_file = os.path.join(args.out_dir,each)
    if os.path.isfile(old_file):
        each_flag,each_name = get_fq_sample_name(each)
        each_first_dir = os.path.join(args.out_dir,each_name)    
        new_file = os.path.join(each_first_dir,'%s_1_%s.fastq.gz' %(each_name,each_flag))
        if not os.path.exists(each_first_dir):
            python_tools.circ_mkdir_unix('"%s"' % each_first_dir)
        prepare_sh.write('mv "%s" "%s"\n' % (old_file,new_file))
prepare_sh.close()

sample_dict = {}
for each_dir in data_dir_list:
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

if os.path.isfile('rawdata_number.json') and os.stat('rawdata_number.json').st_size :
    sample_number_dict = python_tools.load_fn_to_obj('rawdata_number.json')
else :
    sample_number_dict = {}

total_size = []
cp_data_info = open(cp_data_info_file,'w')
for each in sample_dict :
    if each in sample_number_dict :
        sample_dict[each].pre_num = sample_number_dict[each]
    else :
        sample_number_dict[each] = len(sample_dict[each].read1)
    cmd_line = sample_dict[each].add_sup_data(args.out_dir)
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

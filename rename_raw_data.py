'''

rename and zcat rawdata

'''

import os
import sys
import argparse
import python_tools
import RNAseq_tools

parser = argparse.ArgumentParser()
parser.add_argument('--seq_data_dir', help = 'directory store sequencing data', required = True)
parser.add_argument('--seq_dir_name', help = 'sequencing data names, sep with ","', required = True)
parser.add_argument('--analysis_data_dir', help = 'directory store analysis data', required = True)
parser.add_argument('--sample_map', help = 'file map sample id to wgc id', required = True)
args = parser.parse_args()

sample_map_dict = python_tools.table_to_dict(args.sample_map,1,2,False,'\t')

python_tools.circ_mkdir_unix(args.analysis_data_dir)

sample_data_dict = {}
seq_dir_name_list = args.seq_dir_name.split(',')
for each_name in seq_dir_name_list :
    each_dir = os.path.join(args.seq_data_dir,each_name)
    each_dir_files = os.listdir(each_dir) 
    for each_file in each_dir_files :
        each_file_path = os.path.join(each_dir,each_file)
        if os.path.isdir(each_file_path) :
            wgc_id = each_file
            sample_id = sample_map_dict[wgc_id]
            seq_files = os.listdir(each_file_path)
            if sample_id not in sample_data_dict :
                sample_data_dict[sample_id] = RNAseq_tools.rawdata()
                sample_data_dict[sample_id].name = sample_id
            for each_seq_file in seq_files :
                each_seq_path = os.path.join(each_file_path,each_seq_file)
                if each_seq_file.endswith('R1.fastq.gz') :
                    sample_data_dict[sample_id].read1.append(each_seq_path)
                elif each_seq_file.endswith('R2.fastq.gz') :
                    sample_data_dict[sample_id].read2.append(each_seq_path)
                else :
                    pass

merge_cmd = os.path.join(args.analysis_data_dir,'get_analysis_data.sh')
for each in sample_data_dict :
    cmd_line = sample_data_dict[each].merge_rawdata(args.analysis_data_dir)
    python_tools.write_obj_to_file(cmd_line,merge_cmd,True)

os.system('chmod +x %s' % merge_cmd)
python_tools.circ_call_process(merge_cmd)


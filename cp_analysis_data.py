'''
cp analysis data by wgc id

'''

import sys
import os
import argparse
import time
import python_tools
import re

parser = argparse.ArgumentParser()
parser.add_argument('--id_info', help = 'id file', required = True)
parser.add_argument('--data_dir', help = 'data directory', required = True)
parser.add_argument('--out_dir', help = 'table', required = True)
parser.add_argument('--name', help = 'merged table', required = True)
args = parser.parse_args()

id_list = [each.strip() for each in open(args.id_info)]

all_dir_dict = {}
all_ori_file_list =  os.listdir(args.data_dir)
for each in all_ori_file_list :
    each_dir = os.path.join(args.data_dir,each)
    if os.path.isdir(each_dir) :
        wgc_id = re.search(r'(WGC\d{6})',each).groups()[0]
        if wgc_id in id_list :
            all_dir_dict[wgc_id] = each_dir

out_dir = os.path.join(args.out_dir,args.name)
python_tools.circ_mkdir_unix(out_dir)
time_info = time.localtime()
output_time = '%s-%s-%s-%s:%s:%s' % (time_info.tm_year,time_info.tm_mon,time_info.tm_mday,time_info.tm_hour,time_info.tm_min,time_info.tm_sec)
cp_cmd = os.path.join(out_dir,'%s_cp.sh' % output_time)
with open(cp_cmd,'w') as cp_cmd_info :
    for each in all_dir_dict :
        data_dir = all_dir_dict[each]
        cp_cmd_info.write('echo "#### cp %s start ####"\n' % each)
        cp_cmd_info.write('cp -r "%s" %s/%s\n' % (data_dir,out_dir,each))
        cp_cmd_info.write('echo "#### cp %s end ####"\n' % each)
os.system('chmod +x %s' % cp_cmd)
python_tools.circ_call_process(cp_cmd)

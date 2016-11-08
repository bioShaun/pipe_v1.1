'''

get data project info in disk (sequenced in YMKD)

'''

import sys
import os
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--disk_dir',help='data directory ","',required=True)
parser.add_argument('--info_file', help = 'all sample info file', required = True)
parser.add_argument('--disk_info', help = 'output info', required=True)
args = parser.parse_args()

## store sample info
info_dict = {}
with open(args.info_file) as info_file_info :
    for n,eachline in enumerate(info_file_info) :
        if n != 0 :
            eachline_info = eachline.strip().split('\t')
            wgc_id = eachline_info[1]
            proj_id = eachline_info[0]
            info_dict[wgc_id] = proj_id

## acquire sample id
all_file_list = os.listdir(args.disk_dir)

all_dir_list = []
for each in all_file_list :
    file_path = os.path.join(args.disk_dir,each)
    if os.path.isdir(file_path) :
        all_dir_list.append(each)

id_list = []
for each in all_dir_list :
    if 'WGC' in each :
        wgc_id = re.search(r'(WGC\d{6})',each).groups()[0]
        id_list.append(wgc_id)
    else :
        print '%s is not a standard wgc id' % each

## get sample info in disk
with open(args.disk_info,'w') as disk_info_file :
    for each in id_list :
        proj_id = info_dict[each]
        disk_info_file.write('%s\t%s\n' % (each,proj_id))






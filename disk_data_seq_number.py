## coding=utf-8

'''

check sample sequencing number

'''

import sys
import os
import argparse
import re
import json

parser = argparse.ArgumentParser()
parser.add_argument('--sample_dir', help = 'rawdata directory', required = True)
parser.add_argument('--disk_name', help = 'hard disk name', required = True)
args = parser.parse_args()

all_ori_file_list =  os.listdir(args.sample_dir)
sequencing_info_json_file = '/home/lxgui/Documents/label_data/seq_info.json'
output_file = '/home/lxgui/Documents/label_data/seq_info.txt'

if os.path.isfile(sequencing_info_json_file) :
    with open(sequencing_info_json_file) as sequencing_info_json :
        sequencing_info_dict = json.load(sequencing_info_json)
else :
    sequencing_info_dict = {}

sequencing_dict = {}
sequencing_number_file = '/home/lxgui/Documents/paiji/check_up_dowm/0503_all_下机_info_utf8.txt'
with open(sequencing_number_file) as sequencing_number_file_info :
    for n,eachline in enumerate(sequencing_number_file_info) :
        if n != 0 :
            eachline_info = eachline.strip().split('\t')
            wgc_id = eachline_info[0]
            sequencing_number = eachline_info[-1]
            if sequencing_number not in sequencing_dict :
                sequencing_dict[sequencing_number] = [wgc_id]
            else :
                sequencing_dict[sequencing_number].append(wgc_id)
            if wgc_id in sequencing_info_dict and sequencing_number in sequencing_info_dict[wgc_id] :
                pass 
            else :
                sequencing_info_dict.setdefault(wgc_id,{})[sequencing_number] = [eachline,[]]

all_dir_list = []
for each in all_ori_file_list :
    each_dir = os.path.join(args.sample_dir,each)
    if os.path.isdir(each_dir) :
        wgc_id = re.search(r'(WGC\d{6})',each).groups()[0]
        all_dir_list.append(wgc_id)

all_dir_set = set(all_dir_list)
seq_number = None
for each_num in sequencing_dict :
    seq_sample_set = set(sequencing_dict[each_num])
    if all_dir_set.issubset(seq_sample_set) :
        seq_number = each_num
        if all_dir_set.issuperset(seq_sample_set) :
            print 'complete sequencing samples in : %s' % each_num
        else :
            print 'partial sequencing samples in : %s' % each_num

for each_id in all_dir_list :
    seq_num_list = sequencing_info_dict[each_id][seq_number][1]
    if args.disk_name not in seq_num_list :
        sequencing_info_dict[each_id][seq_number][1].append(args.disk_name)

with open(sequencing_info_json_file,'w') as sequencing_info_json :
    json.dump(sequencing_info_dict,sequencing_info_json)

with open(output_file,'w') as output_file_info :
    for each_id in sequencing_info_dict :
        each_seq_num_list = sequencing_info_dict[each_id].keys()
        each_disk_name_list = []
        for each_seq_num in sequencing_info_dict[each_id] :
            disk_name = '--'
            if sequencing_info_dict[each_id][each_seq_num][1] :
                disk_name = ','.join(sequencing_info_dict[each_id][each_seq_num][1])
            each_disk_name_list.append(disk_name)
            # each_other_info = sequencing_info_dict[each_id][each_seq_num][0]
        output_seq_num = ','.join(each_seq_num_list)
        output_disk_name = ','.join(each_disk_name_list)
        output_file_info.write('%s\t%s\t%s\n' % (each_id,output_disk_name,output_seq_num))
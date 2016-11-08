# coding=utf-8

import sys
import os
import math

sample_info_file = sys.argv[1]
dup_info_file = sys.argv[2]
out_dir =  sys.argv[3]

sample_id_dict = {}
data_file = '/home/lxgui/Documents/paiji/check_up_dowm/出库样品信息.txt'

with open(sample_info_file) as sample_info_file_info:
    for eachline in sample_info_file_info:
        eachline_info = eachline.strip().split('\t')
        wgc_id = eachline_info[0]
        sub_sample = eachline_info[1]
        sample_id = eachline_info[1].split('_')[0]
        if sample_id not in sample_id_dict:
            sample_id_dict[sample_id] = [[wgc_id],[sub_sample],[0],[0],0,'no']
        else:
            sample_id_dict[sample_id][0].append(wgc_id)
            sample_id_dict[sample_id][1].append(sub_sample)
            sample_id_dict[sample_id][2].append(0)
            sample_id_dict[sample_id][3].append(0)


with open(dup_info_file) as dup_info_file_info:
    for eachline in dup_info_file_info:
        eachline_info = eachline.strip().split()
        sub_sample_id = eachline_info[0]
        sample_id = sub_sample_id.split('_')[0]
        dup_stat = min(float(eachline_info[1]), float(eachline_info[2]))
        sample_num = sample_id_dict[sample_id][1].index(sub_sample_id)
        sample_id_dict[sample_id][2][sample_num] = dup_stat
        if dup_stat < 80:
            sample_id_dict[sample_id][-1] = 'yes'

header = ''
data_file_dict = {}
with open(data_file) as data_file_info:
    for n,eachline in enumerate(data_file_info):
        eachline_info = eachline.strip().split('\t')
        if n != 0:
            wgc_id = eachline_info[1]
            sample_id = eachline_info[4].split('_')[0]
            sample_datat = float(eachline_info[15])
            if sample_id in sample_id_dict and wgc_id in sample_id_dict[sample_id][0]:
                data_file_dict[wgc_id] = eachline_info
                sample_num = sample_id_dict[sample_id][0].index(wgc_id)
                sample_id_dict[sample_id][3][sample_num] = sample_datat
        else:
            header = '\t'.join(eachline_info[0:18])

output_file1 = os.path.join(out_dir,'sample_detail.txt')
output_file2 = os.path.join(out_dir,'paiji.txt')
output_file1_info = open(output_file1,'w')
output_file2_info = open(output_file2,'w')
output_file2_info.write('%s\n' % header)

paiji_sample_dict = {}
for each_sample in sample_id_dict:
    each_sample_info = sample_id_dict[each_sample]
    wgc_out = '\t'.join(each_sample_info[0])
    sample_out = '\t'.join(each_sample_info[1])
    dup_stat_list = [str(each) for each in each_sample_info[2]]
    sample_dat_list = [str(each) for each in each_sample_info[3]]
    dup_stat_out = '\t'.join(dup_stat_list)
    sample_dat_out = '\t'.join(sample_dat_list)
    need_add = each_sample_info[-1]
    better_dup_sample_num = each_sample_info[2].index(max(each_sample_info[2]))
    add_wgc_id = each_sample_info[0][better_dup_sample_num]
    eff_data = 0
    add_data = 0
    for n,each_data in enumerate(each_sample_info[3]):
        eff_data += each_data*each_sample_info[2][n]/100
    if need_add == 'yes':
        add_data = 100 - eff_data
        if add_data <= 0:
            add_data = 0
            need_add = 'no'
        if need_add == 'yes':
            paiji_sample_dict[add_wgc_id] = [eff_data,add_data]
    output_file1_info.write('{each_sample}\t{sample_out}\t{wgc_out}\t{sample_dat_out}\t{dup_stat_out}\t{eff_data}\t{need_add}\t{add_wgc_id}\t{add_data}\n'.format(**locals()))
output_file1_info.close()

for each_wgc in paiji_sample_dict:
    each_wgc_info = data_file_dict[each_wgc]
    each_eff_data = paiji_sample_dict[each_wgc][0]
    each_add_data = paiji_sample_dict[each_wgc][1]
    each_wgc_info[14] = 100
    each_wgc_info[15] = each_eff_data
    each_wgc_info[16] = each_add_data
    each_wgc_info[17] = math.ceil(each_add_data)
    each_wgc_info_out = [str(each) for each in each_wgc_info]
    each_out = '\t'.join(each_wgc_info_out[0:18])
    output_file2_info.write('%s\n' % each_out)
output_file2_info.close()
        

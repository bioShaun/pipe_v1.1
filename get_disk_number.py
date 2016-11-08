import sys
import os

seq_dict = {}
with open('seq_info.txt') as seq_info_file :
    for eachline in seq_info_file :
        eachline_info = eachline.strip().split('\t')
        disk_list = eachline_info[1].split(',')
        number_list = eachline_info[2].split(',')
        for n,each in enumerate(number_list) :
            if each in seq_dict :
                if disk_list[n] not in seq_dict[each] :
                    seq_dict[each].append(disk_list[n])
            else :
                seq_dict[each] = [disk_list[n]]

sort_number_list = sorted(seq_dict.keys())
with open('seq_num_disk.info','w') as seq_num_disk_info :
    for each in sort_number_list :
        seq_num_disk_info.write('%s\t%s\n' % (each,','.join(seq_dict[each])))
    

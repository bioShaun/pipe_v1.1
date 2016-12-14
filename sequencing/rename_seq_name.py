import sys
import os
import re

if not len(sys.argv) == 3:
    print '    python ' + sys.argv[0] + ' id.map data.dir'
    sys.exit(1)

id_map_file = sys.argv[1]
data_dir = sys.argv[2]
data_type = ('Cleandata', 'Rawdata')

data_dir_list = [os.path.join(data_dir, each) for each in data_type]

id_map_dict = {}
with open(id_map_file) as id_map_file_inf:
    for eachline in id_map_file_inf:
        eachline_inf = eachline.strip().split('\t')
        sample_id = eachline_inf[0]
        lib_id = eachline_inf[1]
        id_map_dict[lib_id] = sample_id

name_map_dict = {}

for each_data_dir in data_dir_list:
    if os.path.isdir(each_data_dir):
        data_file_list = os.listdir(each_data_dir)
        for each_file_dir in data_file_list:
            each_dir_new_name = ""
            for each_id in id_map_dict:
                if each_id in each_file_dir:
                    each_dir_new_name = id_map_dict[each_id]
            each_dir_path = os.path.join(each_data_dir, each_file_dir)
            each_dir_files = os.listdir(each_dir_path)
            for each_file in each_dir_files:
                each_file_new_name = re.sub(each_file_dir, each_dir_new_name, each_file)
                print 'mv {0}/{1} {0}/{2}'.format(each_dir_path, each_file, each_file_new_name)
            print 'mv {0}/{1} {0}/{2}'.format(each_data_dir, each_file_dir, each_dir_new_name)

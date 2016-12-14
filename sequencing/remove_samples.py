import sys
import os

if not len(sys.argv) == 3:
    print '     python ' + sys.argv[0] + ' sample.list data.dir'
    sys.exit(1)

sample_list_file = sys.argv[1]
data_dir = sys.argv[2]
data_type = ('Cleandata', 'Rawdata')

sample_list = [each.strip() for each in open(sample_list_file)]
data_dir_list = [os.path.join(data_dir, each) for each in data_type]

for each_data_dir in data_dir_list:
    if os.path.isdir(each_data_dir):   
        data_file_list = os.listdir(each_data_dir)
        for each_file in data_file_list:
            if each_file in sample_list:
                each_file_path = os.path.join(each_data_dir, each_file)
                print 'rm -r {0}'.format(each_file_path)

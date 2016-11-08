import sys
import os
import random
import math

if not len(sys.argv) == 4:
    print '    python ' + sys.argv[0] + ' data.size sample.size sample.prefix'
    sys.exit(0)

data_size = int(sys.argv[1])
sample_size = int(sys.argv[2])
sample_prefix = sys.argv[3]

def get_data_pool(data_size, data_top = None):
    data_bottom = data_size
    if not data_top:
        data_top = data_size*1.1
    data_list = []
    start_data = data_bottom
    while start_data < data_top:
        start_data += 0.01
        data_list.append(start_data)
    return data_list

def get_sample_names(data_size, sample_prefix):
    data_prefix = int(math.log10(data_size))
    sample_list = []
    for n in range(1, (data_size+1)):
        data_weishu = int(math.log10(n))
        data_name = '{0}{1}{2}'.format(sample_prefix, '0'*(data_prefix-data_weishu), n)
        sample_list.append(data_name)
    return sample_list

data_size_list = random.sample(get_data_pool(data_size), sample_size)
q30_list = random.sample(get_data_pool(90,92), sample_size)
sample_list = get_sample_names(sample_size, sample_prefix)

for n, each in enumerate(sample_list):
    each_data_size = data_size_list[n]
    each_q30 = q30_list[n]
    each_read_num = each_data_size/0.15
    print '{each}\t{each_read_num:.2f}\t{each_data_size:.2f}\tPE150\t{each_q30:.2f}'.format(**locals())


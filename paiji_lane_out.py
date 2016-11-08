#! /usr/bin/python
#coding=utf-8

import sys
import os
reload(sys)
sys.setdefaultencoding('utf8')

if not len(sys.argv) == 4:
    print '    python ' + sys.argv[0] + ' lane.inf seq.data out.dir'
    sys.exit(0)

lane_inf_file = sys.argv[1]
seq_data = sys.argv[2]
out_dir = sys.argv[3]

lane_dict = {}
lane_inf_data_file = os.path.join(out_dir, 'lane.data.txt')
lane_summary = os.path.join(out_dir, 'lane.summary.txt')

seq_data_dict = {}
with open(seq_data) as seq_data_inf:
    for n, eachline in enumerate(seq_data_inf):
        if n != 0 :
            eachline_inf = eachline.strip().split('\t')
            wgc_id = eachline_inf[0].strip()
            data_size = eachline_inf[-2].strip()
            seq_num = eachline_inf[-1].strip()
            seq_data_dict.setdefault(wgc_id, {})[seq_num] = data_size

lane_inf_data_file_inf =  open(lane_inf_data_file, 'w')
with open(lane_inf_file) as lane_inf:
    for n, eachline in enumerate(lane_inf):
        eachline = eachline.strip()
        if n == 0:
            lane_inf_data_file_inf.write('%s\tData_size\n' % eachline)
        else:
            eachline_inf = eachline.split('\t')
            wgc_id = eachline_inf[0].strip()
            seq_num = eachline_inf[-1].strip()
            lane_num = eachline_inf[1].strip()
            data_size = float(seq_data_dict[wgc_id][seq_num])
            lane_inf_data_file_inf.write('%s\t%s\n' % (eachline, data_size))
            if seq_num in lane_dict and lane_num in lane_dict[seq_num]:
                lane_dict[seq_num][lane_num] += data_size
            else:
                lane_dict.setdefault(seq_num, {})[lane_num] = data_size

lane_data_120_list = []
lane_data_110_list = []
lane_data_lt_110_list = []

for each_seq in lane_dict:
    for each_lane in lane_dict[each_seq]:
        each_data_size = lane_dict[each_seq][each_lane]
        each_lane_detail = '%s:lane%s(%sG)' % (each_seq, each_lane, each_data_size)
        if each_data_size >= 120:
            lane_data_120_list.append(each_lane_detail)
        elif each_data_size >= 110:
            lane_data_110_list.append(each_lane_detail)
        else:
            lane_data_lt_110_list.append(each_lane_detail)


with open(lane_summary, 'w') as lane_summary_inf:
    lane_summary_inf.write('Lane数据量大于120G数量：%s\n' % len(lane_data_120_list))
    lane_summary_inf.write('Lane数据量小于120G，大于110G数量：%s\n' % len(lane_data_110_list))
    lane_summary_inf.write('Lane数据量小于110G数量：%s\n' % len(lane_data_lt_110_list))
    lane_summary_inf.write('Lane数据量大于120G详情：%s\n' % ','.join(lane_data_120_list))
    lane_summary_inf.write('Lane数据量小于120G，大于110G详情：%s\n' % ','.join(lane_data_110_list))
    lane_summary_inf.write('Lane数据量小于110G详情：%s\n' % ','.join(lane_data_lt_110_list))


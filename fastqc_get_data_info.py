import sys
import os

if not len(sys.argv) == 3:
    print 'python ' + sys.argv[0] + ' sample_list qc_dir'
    sys.exit(0)

sample_list_file = sys.argv[1]
qc_dir = sys.argv[2]

sample_list = [each.strip() for each in open(sample_list_file)]
sample_info_dict = {}

for each_sample in sample_list:
    sample_info_dict[each_sample] = [0, 0, 0, []]
    for n in (1, 2):
        each_qc_dir = os.path.join(qc_dir, '%s_%s.fq_fastqc' % (each_sample, n))
        each_qc_file = os.path.join(each_qc_dir, 'fastqc_data.txt')
        with open(each_qc_file) as each_qc_file_info:
            q30_flag = 0
            for eachline in each_qc_file_info:
                eachline_info = eachline.strip().split('\t')
                if 'Sequence length' in eachline:                    
                    sample_info_dict[each_sample][0] = int(eachline_info[1])
                if 'Total Sequences' in eachline:
                    sample_info_dict[each_sample][1] += int(eachline_info[1])
                if 'Per sequence quality scores' in eachline:
                    q30_flag = 1
                    continue
                if eachline.startswith('>>'):
                    q30_flag = 0
                if not eachline.startswith('#') and q30_flag == 1 and int(float(eachline_info[0])) >= 30 :
                    sample_info_dict[each_sample][2] += int(float(eachline_info[1]))
                if '#Total Duplicate Percentage' in eachline:
                    dup = round(float(eachline_info[1]), 2)
                    sample_info_dict[each_sample][3].append('%s(read%s)' % (dup, n))

qc_results = os.path.join(qc_dir, 'qc_results.txt')
qc_results_info = open(qc_results, 'w')
qc_results_info.write('Sample_ID\tReads_number(M)\tData_size(G)\tQ30(%)\tDuplication(%)\n')
for each_sample in sample_list:
    reads_num = round(sample_info_dict[each_sample][1]/float(10**6), 2)
    data_size = round(sample_info_dict[each_sample][1]*sample_info_dict[each_sample][0]/float(1000**3), 2)
    q30 = round(sample_info_dict[each_sample][2]*100/float(sample_info_dict[each_sample][1]), 2)
    dup = ','.join(sample_info_dict[each_sample][3])
    qc_results_info.write('%s\t%s\t%s\t%s\t%s\n' % (each_sample, reads_num, data_size, q30, dup))
qc_results_info.close()






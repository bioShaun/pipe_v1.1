import sys
import os

if not len(sys.argv) == 3:
    print 'python ' + sys.argv[0] + ' sample_list qc_dir > out.summary'
    sys.exit(0)

sample_list_file = sys.argv[1]
qc_dir = sys.argv[2]

sample_list = [each.strip() for each in open(sample_list_file)]
sample_info_dict = {}

for each_sample in sample_list:
    sample_info_dict[each_sample] = [0, 0, 0, []]
    for n in (1, 2):
        each_qc_dir = os.path.join(qc_dir, '%s_%s_clean.fq_fastqc' % (each_sample, n))
        each_qc_file = os.path.join(each_qc_dir, 'fastqc_data.txt')
        with open(each_qc_file) as each_qc_file_info:
            q30_flag = 0
            for eachline in each_qc_file_info:
                eachline_info = eachline.strip().split('\t')
                if 'Sequence length' in eachline:                    
                    sample_info_dict[each_sample][0] = int(eachline_info[1])
                if 'Total Sequences' in eachline:
                    sample_info_dict[each_sample][1] += int(eachline_info[1])
                if eachline.startswith("%GC"):
                    gc_content = int(eachline_info[1])
                    sample_info_dict[each_sample][3].append(gc_content)                   
                if 'Per sequence quality scores' in eachline:
                    q30_flag = 1
                    continue
                if eachline.startswith('>>'):
                    q30_flag = 0
                if not eachline.startswith('#') and q30_flag == 1 and int(float(eachline_info[0])) > 30 :
                    sample_info_dict[each_sample][2] += int(float(eachline_info[1]))

print 'Sample_ID\tReads_number(M)\tReads_length(bp)\tData_size(G)\tQ30(%)\tGC(%)'
for each_sample in sample_list:
    reads_num = round(sample_info_dict[each_sample][1]/float(10**6), 2)
    read_length = sample_info_dict[each_sample][0] 
    data_size = round(sample_info_dict[each_sample][1]*sample_info_dict[each_sample][0]/float(1000**3), 2)
    q30 = round(sample_info_dict[each_sample][2]*100/float(sample_info_dict[each_sample][1]), 2)
    gc_content = sum(sample_info_dict[each_sample][3])/float(len(sample_info_dict[each_sample][3]))
    print '%s\t%s\t%s\t%s\t%s\t%s' % (each_sample, reads_num, read_length, data_size, q30, gc_content)






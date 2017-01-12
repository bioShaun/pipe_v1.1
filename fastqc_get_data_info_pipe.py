import sys
import os
import python_tools

if not len(sys.argv) == 4:
    print 'python ' + sys.argv[0] + ' sample_list qc_dir out_dir'
    sys.exit(0)

sample_list_file = sys.argv[1]
qc_dir = sys.argv[2]
out_dir = sys.argv[3]

sample_list = [each.strip() for each in open(sample_list_file)]
sample_info_dict = {}
reads_quality_dir = os.path.join(out_dir, 'reads_quality')
python_tools.circ_mkdir_unix(reads_quality_dir)
merged_quality_file = os.path.join(reads_quality_dir, 'all.reads_quality.txt')
merged_quality_file_inf = open(merged_quality_file, 'w')
merged_quality_file_inf.write('Sample_ID\tQuality\tCount\tPercent\n')

for each_sample in sample_list:
    sample_info_dict[each_sample] = [0, 0, 0, []]
    reads_quality_dict = {}
    each_reads_quality_file = os.path.join(reads_quality_dir, '{}_reads_quality.txt'.format(each_sample))
    for n in (1, 2):
        each_qc_dir = os.path.join(qc_dir, '%s_%s.clean.fq_fastqc' % (each_sample, n))
        each_qc_file = os.path.join(each_qc_dir, 'fastqc_data.txt')
        with open(each_qc_file) as each_qc_file_info:
            q30_flag = 0
            for eachline in each_qc_file_info:
                eachline_info = eachline.strip().split('\t')
                if 'Sequence length' in eachline:                    
                    sample_info_dict[each_sample][0] = int(eachline_info[1])
                if 'Total Sequences' in eachline:
                    sample_info_dict[each_sample][1] += int(eachline_info[1])
                if eachline.startswith('%GC'):
                    gc_content = int(eachline_info[1])
                    sample_info_dict[each_sample][3].append(gc_content)                   
                if 'Per sequence quality scores' in eachline:
                    q30_flag = 1
                    continue
                if eachline.startswith('>>'):
                    q30_flag = 0
                    quality_flag = 0                    
                if not eachline.startswith('#') and q30_flag == 1:
                    each_quality = int(eachline_info[0])
                    each_count = float(eachline_info[1])
                    if each_sample in reads_quality_dict and each_quality in reads_quality_dict[each_sample]:
                        reads_quality_dict[each_sample][each_quality] += each_count
                    else:
                        reads_quality_dict.setdefault(each_sample, {})[each_quality] = each_count
                    if int(float(eachline_info[0])) > 30 :
                        sample_info_dict[each_sample][2] += int(float(eachline_info[1]))
    with open(each_reads_quality_file, 'w') as each_reads_quality_file_inf:
        quality_list = sorted(reads_quality_dict[each_sample].keys())
        each_reads_quality_file_inf.write('Quality\tCount\tProportion\n')
        quality_count_list = []
        for each_quality in quality_list:
            quality_count_list.append(reads_quality_dict[each_sample][each_quality])
        for n, each_quality in enumerate(quality_list):
            each_quality_portion = float(quality_count_list[n])/sum(quality_count_list)
            each_quality_count = quality_count_list[n]
            each_reads_quality_file_inf.write('{each_quality}\t{each_quality_count}\t{each_quality_portion}\n'.format(**locals()))
            merged_quality_file_inf.write('{each_sample}\t{each_quality}\t{each_quality_count}\t{each_quality_portion}\n'.format(**locals()))
merged_quality_file_inf.close()

qc_summary_file = os.path.join(out_dir, 'qc.summary.txt')
with open(qc_summary_file, 'w') as qc_summary_file_inf:
    qc_summary_file_inf.write('Sample_ID\tReads_number(M)\tReads_length(bp)\tData_size(G)\tQ30(%)\tGC(%)\n')
    for each_sample in sample_list:
        reads_num = round(sample_info_dict[each_sample][1]/float(10**6), 2)
        read_length = sample_info_dict[each_sample][0] 
        data_size = round(sample_info_dict[each_sample][1]*sample_info_dict[each_sample][0]/float(1000**3), 2)
        q30 = round(sample_info_dict[each_sample][2]*100/float(sample_info_dict[each_sample][1]), 2)
        gc_content = sum(sample_info_dict[each_sample][3])/float(len(sample_info_dict[each_sample][3]))
        qc_summary_file_inf.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (each_sample, reads_num, read_length, data_size, q30, gc_content))






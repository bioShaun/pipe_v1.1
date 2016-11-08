from __future__ import division
import sys
import os
import python_tools

if not len(sys.argv) == 4:
    print '    print ' + sys.argv[0] + ' sRNA.analsys.dir sRNA.qc.summary sRNA.length.dir'
    sys.exit(1)

sRNA_analysis_dir = sys.argv[1]
out_qc_summary = sys.argv[2]
sRNA_length_dir = sys.argv[3]

python_tools.circ_mkdir_unix(sRNA_length_dir)
sRNA_ori_summary_file = os.path.join(sRNA_analysis_dir,'SampleSummary.xls')
sRNA_data_dict = {}

with open(sRNA_ori_summary_file) as sRNA_ori_summary_file_inf:
    for n, eachline in enumerate(sRNA_ori_summary_file_inf):
        if n != 0:
            eachline_inf = eachline.strip().split('\t')
            sample_id = eachline_inf[0]
            total_reads = int(eachline_inf[4])
            mapped_reads = int(eachline_inf[5])
            # mapping_rate = 100*round(mapped_reads/mapped_reads, 2)
            # mapping_rate_rep = '%s%%' mapping_rate
            sRNA_data_dict[sample_id] = [total_reads, mapped_reads]

sample_length_merged_file = os.path.join(sRNA_length_dir, 'Sample.length.txt')
sample_length_merged_file_inf = open(sample_length_merged_file, 'w')
sample_length_merged_file_inf.write('Sample\tLength\tCount\n')
for each_sample in sRNA_data_dict:
    each_sample_length_file = os.path.join(sRNA_length_dir,'%s.length.txt' % each_sample)
    each_sample_length_file_inf = open(each_sample_length_file, 'w')
    each_sample_qc_file = os.path.join(sRNA_analysis_dir, 'qc/fastqc_posttrim/%s/%s.cutadapt_fastqc/fastqc_data.txt' % (each_sample, each_sample))
    with open(each_sample_qc_file) as each_sample_qc_file_inf:
        quality_flag = 0
        length_flag = 0
        q30_reads = 0
        for eachline in each_sample_qc_file_inf:
            if eachline.startswith('#Quality'):
                quality_flag = 1
                continue
            if eachline.startswith('#Length'):
                length_flag = 1
                eachline = eachline.lstrip('#')
                each_sample_length_file_inf.write(eachline)
                continue
            if eachline.startswith('>>END_MODULE'):
                quality_flag = 0
                length_flag = 0
            if length_flag:
                sample_length_merged_file_inf.write('{each_sample}\t{eachline}'.format(**locals()))
                each_sample_length_file_inf.write(eachline)
            if quality_flag:
                eachline_inf = eachline.strip().split('\t')
                quality_score = int(eachline_inf[0])
                readcount = float(eachline_inf[1])
                if quality_score > 30:
                    q30_reads += readcount
        sRNA_data_dict[each_sample].append(q30_reads)
        each_sample_length_file_inf.close()

with open(out_qc_summary, 'w') as out_qc_summary_inf:
    sample_list = sorted(sRNA_data_dict.keys())
    out_qc_summary_inf.write('Sample_ID\tClean_reads\tMapped_reads\tMapping_rate\tQ30\n')    
    for each_sample in sample_list:
        total_reads, mapped_reads, q30_reads = sRNA_data_dict[each_sample]
        mapping_rate = 100*mapped_reads/total_reads
        q30_rate = 100*q30_reads/total_reads
        out_qc_summary_inf.write('{each_sample}\t{total_reads:,}\t{mapped_reads:,}\t{mapping_rate:.2f}%\t{q30_rate:.2f}%\n'.format(**locals()))




import os
import sys

if not len(sys.argv) == 3:
    print 'python ' + sys.argv[0] + ' sample_list quant_dir'
    sys.exit(0)

sample_list_file = sys.argv[1]
quant_dir = sys.argv[2]

sample_list = [each.strip() for each in open(sample_list_file)]
quant_dict = {}

for each_sample in sample_list:
    each_sample_dir = os.path.join(quant_dir, each_sample)
    each_sample_quant = os.path.join(each_sample_dir, 'miRNAs_expressed_all_samples_now.csv')
    each_sample_quant_info = open(each_sample_quant)
    for n, eachline in enumerate(each_sample_quant_info):
        if n != 0:
            eachline_info = eachline.strip().split('\t')
            miRNA_id = eachline_info[0]
            count = int(float(eachline_info[1]))
            if miRNA_id in quant_dict and each_sample in quant_dict[miRNA_id]:
                quant_dict[miRNA_id][each_sample] += count
            else:
                quant_dict.setdefault(miRNA_id, {})[each_sample] = count

merged_table = os.path.join(quant_dir, 'miRNA_expression_count.txt')
with open(merged_table, 'w') as merged_table_info:
    merged_table_info.write('miRNA\t%s\n' % '\t'.join(sample_list))
    for each_miRNA in quant_dict:
        quant_list = []
        for each_sample in sample_list:
            quant_list.append(str(quant_dict[each_miRNA][each_sample]))       
        merged_table_info.write('%s\t%s\n' % (each_miRNA, '\t'.join(quant_list)))


'''
merge novel and known miRNA count table

'''

import sys
import os

if not len(sys.argv) == 4:
    print 'python ' + sys.argv[0] + ' mature_miRNA_expression.xls novel_miRNA.xls merged.table.xls'
    sys.exit(0)

known_miRNA = sys.argv[1]
novel_miRNA = sys.argv[2]
merge_table = sys.argv[3]

novel_miRNA_dir, novel_miRNA_name = os.path.split(novel_miRNA)
new_novel_miRNA = os.path.join(novel_miRNA_dir,'novel_miRNA_detail.txt')

novel_miRNA_list = []
novel_miRNA_dict = {}
sample_index_dict = {}

new_novel_miRNA_info = open(new_novel_miRNA, 'w')
with open(novel_miRNA) as novel_miRNA_info:
    for n, eachline in enumerate(novel_miRNA_info):
        eachline_info = eachline.strip().split('\t')
        if n == 0:
            new_novel_miRNA_info.write('ID\t%s' % eachline)
            for m, each_col in enumerate(eachline_info):
                if 'ReadCount' in each_col:
                    sample_id = each_col.split('_ReadCount')[0]
                    sample_index_dict[m] = sample_id
        else:            
            novel_miRNA_id = 'novel_%s' % n
            novel_miRNA_list.append(novel_miRNA_id)
            new_novel_miRNA_info.write('%s\t%s' % (novel_miRNA_id, eachline))
            for m, each_col in enumerate(eachline_info):
                if m in sample_index_dict:
                    sample_id = sample_index_dict[m]
                    if each_col == '-':
                        exp_value = "0"
                    else:
                        each_col_list = each_col.split(',')
                        each_col_list = [int(each.strip()) for each in each_col_list]
                        exp_value = str(sum(each_col_list))
                    novel_miRNA_dict.setdefault(novel_miRNA_id, {})[sample_id] = exp_value

sample_list = []
merge_table_info = open(merge_table, 'w')
with open(known_miRNA) as known_miRNA_info:
    for n, eachline in enumerate(known_miRNA_info):
        eachline_info = eachline.strip().split('\t')
        if n == 0:
            merge_table_info.write(eachline)
            eachline_info.pop(0)
            sample_list = eachline_info
        else:
            merge_table_info.write(eachline)

for each_novel in novel_miRNA_list:
    each_novel_info = [each_novel]
    for each_sample in sample_list:
        each_novel_info.append(novel_miRNA_dict[each_novel][each_sample])
    each_novel_out = '\t'.join(each_novel_info)
    merge_table_info.write('%s\n' % each_novel_out)
merge_table_info.close()



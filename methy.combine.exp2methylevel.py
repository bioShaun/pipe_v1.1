import sys
import os

if not len(sys.argv) == 4:
    print '    python ' + sys.argv + ' exp.table meth.dir merged.table'
    sys.exit(1)

exp_table = sys.argv[1]
meth_dir = sys.argv[2]
merged_table = sys.argv[3]

gene_region = ('promoter', 'genebody', 'tts')

sample_inf_dict = {}

sample_list = []
with open(exp_table) as exp_table_inf:
    for n, eachline in enumerate(exp_table_inf):
        eachline_inf = eachline.strip().split('\t')
        if n == 0:
            sample_list = eachline_inf[1:]
        else:
            gene_id = eachline_inf[0]
            for m, each_exp in enumerate(eachline_inf[1:]):
                sample_name = sample_list[m]
                sample_inf_dict.setdefault(gene_id, {})[sample_name] = ['0', '0', '0', each_exp]

for each_sample in sample_list:
    each_methy_file = os.path.join(meth_dir, '%s.gene.meth.level' % each_sample)
    with open(each_methy_file) as each_methy_file_inf:
        for eachline in each_methy_file_inf:
            eachline_inf = eachline.strip().split('\t')
            gene_id, region = eachline_inf[3].split('|')
            methy_level = eachline_inf[-1]
            region_index = gene_region.index(region)
            if gene_id in sample_inf_dict and methy_level != '.' :
                sample_inf_dict[gene_id][each_sample][region_index] = methy_level

gene_feature = list(gene_region)
gene_feature.append('expression')
table_header = ['%s_%s' % (a, b) for b in gene_feature for a in sample_list]
table_header.insert(0, 'Gene_ID')

with open(merged_table, 'w') as merged_table_inf:
    merged_table_inf.write('%s\n' % '\t'.join(table_header))
    for each_gene in sample_inf_dict:
        out_list = []
        for each_sample in sample_list:
            out_list.append(sample_inf_dict[each_gene][each_sample])
        merged_list = []
        for i in range(len(gene_feature)):
            for j in range(len(sample_list)):
                merged_list.append(out_list[j][i])
        merged_table_inf.write('%s\n' % '\t'.join(merged_list))
        



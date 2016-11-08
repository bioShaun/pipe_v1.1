import sys
import os
import python_tools
from HTSeq import GFF_Reader

if not len(sys.argv) == 3:
    print "    python " + sys.argv[0] + ' filter.inf filter.stat'
    sys.exit(0)

def get_gene_number(gtf):
    gene_dict = {}
    for eachline in GFF_Reader(gtf):
        gene_id = eachline.attr['gene_id']
        if gene_id not in gene_dict:
            gene_dict[gene_id] = 1
    return len(gene_dict.keys())

filter_inf = sys.argv[1]
filter_stat = sys.argv[2]

filter_gtf_dict = python_tools.table_to_dict(filter_inf,1,2,False)
gene_num_dict = {}
for each in filter_gtf_dict:
    each_gene_number = get_gene_number(filter_gtf_dict[each])
    gene_num_dict[each] = each_gene_number

sorted_gene_num = sorted(gene_num_dict.iteritems(), key=lambda gene_num_dict : gene_num_dict[1])

with open(filter_stat, 'w') as filter_stat_inf:
    for each in sorted_gene_num:
        filter_stat_inf.write('%s\t%s\n' % (each[0], each[1]))


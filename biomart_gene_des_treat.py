import sys
import os
import csv

if not len(sys.argv) == 3:
    print '    python ' + sys.argv[0] + ' biomart.download treat.biomart.anno'
    sys.exit(1)

biomart_dl = sys.argv[1]
treat_out = sys.argv[2]

reader = csv.reader(file(biomart_dl, 'rb'))
treat_out_inf = open(treat_out, 'w')

gene_anno_dict = {}    
for n, each_record in enumerate(reader):
    if n == 0:
        header = '\t'.join(each_record)
        treat_out_inf.write('{header}\n'.format(**locals()))
    else:
        gene_id = each_record[0]
        gene_inf = each_record[1:]
        if gene_id not in gene_anno_dict:
            gene_anno_dict[gene_id] = []
            for i in range(len(gene_inf)):
                gene_anno_dict[gene_id].append([])
        for m, each_inf in enumerate(gene_inf):
            if each_inf and each_inf not in gene_anno_dict[gene_id][m]:
                gene_anno_dict[gene_id][m].append(each_inf)

for each_gene in gene_anno_dict:
    out_list = []
    for each_inf in gene_anno_dict[each_gene]:
        if each_inf:
            out_list.append('||'.join(each_inf))
        else:
            out_list.append('--')
    out_inf = '\t'.join(out_list)
    treat_out_inf.write('{each_gene}\t{out_inf}\n'.format(**locals()))
treat_out_inf.close()

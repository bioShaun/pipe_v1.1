import sys
import os
import RNAseq_tools

if not len(sys.argv) == 3:
    print '    python ' + sys.argv[0] + ' kobas.kegg.annotate kegg.annotate.txt'
    sys.exit(1)

kobas_anno = sys.argv[1]
extract_out = sys.argv[2]

gene_ko_dict = {}

with open(kobas_anno) as kobas_anno_inf:
    gene_id = ""
    for eachline in kobas_anno_inf:
        eachline_inf = eachline.rstrip().split('\t')
        if '#' not in eachline_inf[0] and 'Query' in eachline_inf[0]:
            gene_id = eachline_inf[1]
        if len(eachline_inf) == 4 and eachline_inf[2] == 'KEGG PATHWAY':
            ko_id = eachline_inf[-1]
            ko_des = eachline_inf[1]
            if gene_id not in gene_ko_dict:
                gene_ko_dict[gene_id] = [[ko_id], [ko_des]]
            else:
                if ko_id not in gene_ko_dict[gene_id][0]:
                    gene_ko_dict[gene_id][0].append(ko_id)
                    gene_ko_dict[gene_id][1].append(ko_des)

with open(extract_out, 'w') as extract_out_inf:
    extract_out_inf.write('Gene_ID\tKO_ID\tKO_description\n')
    for each_gene in gene_ko_dict:
        ko_id = '||'.join(gene_ko_dict[each_gene][0])
        ko_des = '||'.join(gene_ko_dict[each_gene][1])
        extract_out_inf.write('{each_gene}\t{ko_id}\t{ko_des}\n'.format(**locals()))




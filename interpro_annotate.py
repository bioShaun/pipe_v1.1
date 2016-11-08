import sys
import os
import python_tools

if not len(sys.argv) == 4:
    print '    python ' + sys.argv[0] + ' exp.table gene.interpro.annotate exp.annotate.table'
    sys.exit(0)

exp_table = sys.argv[1]
interpro_json = sys.argv[2]
exp_anno = sys.argv[3]

interpro_dict = python_tools.load_fn_to_obj(interpro_json)
exp_anno_inf = open(exp_anno, 'w')
with open(exp_table) as exp_table_inf:
    for n, eachline in enumerate(exp_table_inf):
        eachline = eachline.strip()
        if n == 0:
            exp_anno_inf.write('%s\tGene_name\tInterpro_ID\tInterpro_description\n' % eachline)
        else:
            eachline_inf = eachline.split('\t')
            gene_id = eachline_inf[0]
            annotate_list = ['--']*3
            if gene_id in interpro_dict:
                for m, each in enumerate(annotate_list):
                    if interpro_dict[gene_id][m]:
                        annotate_list[m] = ','.join(interpro_dict[gene_id][m])
            annotate_out = '\t'.join(annotate_list)
            exp_anno_inf.write('%s\t%s\n' % (eachline, annotate_out))



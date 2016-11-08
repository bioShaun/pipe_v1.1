import sys
import os
from HTSeq import GFF_Reader
import python_tools

if not len(sys.argv) == 4:
    print '    python ' + sys.argv[0] + ' lncRNA.feature novel.gtf add.novel.gtf'
    sys.exit(0)

lncRNA_feature = sys.argv[1]
novel_gtf = sys.argv[2]
add_gtf = sys.argv[3]

lncRNA_tr_dict = {}
lncRNA_gene_dict = {}
with open(lncRNA_feature) as lncRNA_feature_inf:
    for n, eachline in enumerate(lncRNA_feature_inf):
        if n != 0:
            eachline_inf = eachline.strip().split('\t')
            tr_id = eachline_inf[4]
            gene_id = eachline_inf[5]
            tr_type = eachline_inf[-1]
            lncRNA_tr_dict[tr_id] = tr_type
            lncRNA_gene_dict[gene_id] = tr_type

out_list = []
for eachline in GFF_Reader(novel_gtf):
    gene_id = eachline.attr['gene_id']
    transcript_id  = eachline.attr['transcript_id']
    gene_type = tr_type = 'TUCP'
    if gene_id in lncRNA_gene_dict:
        gene_type = lncRNA_gene_dict[gene_id]
    if transcript_id in lncRNA_tr_dict:
        tr_type = lncRNA_tr_dict[transcript_id]
    out_list.append('%s; gene_biotype "%s"; transcript_biotype "%s";' % (eachline.get_gff_line().strip(), gene_type, tr_type))
python_tools.write_obj_to_file(out_list, add_gtf)






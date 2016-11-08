import sys
import os
from HTSeq import GFF_Reader

if not len(sys.argv) == 5:
    print '    python ' + sys.argv[0] + ' lncRNA_classify_file miRNA_gene_file input_gtf output_gtf'
    sys.exit(0)

lncRNA_classify_file = sys.argv[1]
miRNA_gene_file = sys.argv[2]
input_gtf = sys.argv[3]
output_gtf = sys.argv[4]

miRNA_gene_list = [each.strip() for each in open(miRNA_gene_file)]

sense_overlapping_tr_dict = {}
with open(lncRNA_classify_file) as lncRNA_classify_file_info:
    for n, eachline in enumerate(lncRNA_classify_file_info):
        if n != 0:
            eachline_info = eachline.strip().split('\t')
            each_tr_id = eachline_info[2]            
            if eachline_info[5] == 'sense' and eachline_info[6] == 'genic' and eachline_info[9] == 'exonic':
                sense_overlapping_tr_dict[each_tr_id] = 'so'
                if eachline_info[3] in miRNA_gene_list:
                    sense_overlapping_tr_dict[each_tr_id] = 'mi'

output_gtf_info = open(output_gtf, 'w')
for eachline in GFF_Reader(input_gtf):
    tr_id = eachline.attr['transcript_id']
    output_line = '%s;' % eachline.get_gff_line().strip()
    if tr_id in sense_overlapping_tr_dict:
        if sense_overlapping_tr_dict[tr_id] == 'mi':
            output_line = '%s transcript_type "miRNA_host";' % output_line
        else:
            continue
    output_gtf_info.write('%s\n' % output_line)
output_gtf_info.close()




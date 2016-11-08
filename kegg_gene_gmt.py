import sys
import os
import argparse
import RNAseq_tools

parser = argparse.ArgumentParser()
parser.add_argument('--kegg_annotate', help = 'go download from biomart.', required = True)
parser.add_argument('--gtf', help = 'gtf file.', required = True)
#parser.add_argument('--species', help = 'maf species.', required = True)
parser.add_argument('--output', help = 'output file', required = True)
args = parser.parse_args()

transcript_dict = RNAseq_tools.get_transcript_info(args.gtf)

kegg_gene_dict = {}
output_info = open(args.output,'w')
with open(args.kegg_annotate) as kegg_annotate_info:
    for eachline in kegg_annotate_info:
        eachline_info = eachline.rstrip().split('\t')
        if 'Query:' in eachline:
            transcript_id = eachline_info[1]
            gene_id = transcript_dict[transcript_id]['gene_id']
        if len(eachline_info) == 4 and eachline_info[2] == 'KEGG PATHWAY':
            pathway_id = eachline_info[-1]
            pathway_name = eachline_info[1]
            if pathway_id not in kegg_gene_dict:
                kegg_gene_dict[pathway_id] = [[gene_id],pathway_name]
            else:
                if gene_id not in kegg_gene_dict[pathway_id][0]:
                    kegg_gene_dict[pathway_id][0].append(gene_id)

for each_id in kegg_gene_dict:
    gene_id_info = '\t'.join(kegg_gene_dict[each_id][0])
    each_name = kegg_gene_dict[each_id][1]
    output_info.write('%s\t%s\t%s\n' % (each_id,each_name,gene_id_info))
output_info.close()

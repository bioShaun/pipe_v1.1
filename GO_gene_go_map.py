import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--biomart_go', help = 'go download from biomart.', required = True)
#parser.add_argument('--species', help = 'maf species.', required = True)
parser.add_argument('--out_dir', help = 'output directory', required = True)
args = parser.parse_args()

gene_go_dict = {}
go_gene_dict = {}

out_dir = os.path.abspath(args.out_dir)

go_file_name = os.path.basename(args.biomart_go)
go_prefix = os.path.splitext(go_file_name)[0]

with open(args.biomart_go) as go_data_info :
    for n,eachline in enumerate(go_data_info) :
        if n != 0 :
            eachline_info = eachline.strip().split(',')
            gene_id = eachline_info[0]
            go_id = eachline_info[1]
            if not go_id:
                continue
            if gene_id not in gene_go_dict :
                gene_go_dict[gene_id] = [go_id]
            else :
                gene_go_dict[gene_id].append(go_id)
            if go_id not in go_gene_dict:
                go_gene_dict[go_id] = [gene_id]
            else:
                go_gene_dict[go_id].append(gene_id)

topgo_go_file = os.path.join(out_dir,'%s_gene_go.txt' % go_prefix)
page_go_file = os.path.join(out_dir,'%s_go.gmt' % go_prefix)
topgo_go_file_info = open(topgo_go_file,'w')
page_go_file_info = open(page_go_file,'w')

for each_gene in gene_go_dict :
    each_go_info = ','.join(gene_go_dict[each_gene])
    topgo_go_file_info.write('%s\t%s\n' % (each_gene,each_go_info))
topgo_go_file_info.close()

for each_go in go_gene_dict:
    each_gene_info = '\t'.join(go_gene_dict[each_go])
    page_go_file_info.write('%s\t%s\t%s\n' % (each_go,each_go,each_gene_info))
page_go_file_info.close()

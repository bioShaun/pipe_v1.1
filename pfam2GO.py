## coding=utf-8

import sys
import os
import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument('--gene_trans',help = 'gene transcript mapping file.', required = True)
parser.add_argument('--pfam_annotate',help = 'pfam annotation file.', required =True)
parser.add_argument('--outgo' , help = 'output go annotation.', default = '1')
args = parser.parse_args()

tr_gene_dict = {}
with open(args.gene_trans) as gene_trans_info:
    for eachline in gene_trans_info:
        eachline_info = eachline.strip().split('\t')
        gene_id = eachline_info[0]
        tr_id = eachline_info[1]
        tr_gene_dict[tr_id] = gene_id

pfam2go_json_file = '/home/public/database/Pfam/pfam2go.json'
with open(pfam2go_json_file) as pfam2go_json_file_info:
    pfam2go_dict = json.load(pfam2go_json_file_info)

outgo_info = open(args.outgo,'w')
outgo_info.write('GeneID,GO\n')

gene_go_dict = {}
with open(args.pfam_annotate) as pfam_annotate_info:
    for eachline in pfam_annotate_info:
        if not eachline.startswith('#') and eachline.strip() != "" :
            each_tr_pfam = eachline.strip().split()
            tr_id = ''.join(each_tr_pfam[0].split('.')[:-1])
            gene_id = tr_gene_dict[tr_id]
            pfam_id = each_tr_pfam[5].split('.')[0]
            if pfam_id in pfam2go_dict:
                go_id_list = pfam2go_dict[pfam_id]
                for each_go in go_id_list:
                   ## 去除冗余结果
                   if gene_id in gene_go_dict and each_go in gene_go_dict[gene_id]:
                       continue
                   else:
                       gene_go_dict.setdefault(gene_id, {})[each_go] = 1
                       outgo_info.write('%s,%s\n' %(gene_id,each_go))                    
outgo_info.close()

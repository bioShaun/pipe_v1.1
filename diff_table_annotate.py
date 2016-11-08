'''
annotate diff table

input1: csv format ensembl biomart annotation ('Gene stable ID', 'GO term accession', 'GO term name', 'InterPro ID', 'InterPro description').
input2: KOBAS annotation
input3: gtf file
input4: diff analysis output

'''

import csv
import json
import argparse
import RNAseq_tools
from goatools import obo_parser
import sys
import python_tools

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
group.add_argument('--biomart', help = 'biomart annotation')
group.add_argument('--go', help = 'go annotation')
parser.add_argument('--kobas', help = 'KOBAS annotation', required = True)
parser.add_argument('--gtf', help = 'gtf format file', required = True)
parser.add_argument('--diff_out', help = 'differential analysis output.', required = True)
parser.add_argument('--output', help = 'output directory', required = True)
parser.add_argument('--type', help = 'analysis type', choices = ['transcript', 'gene'],  default = 'gene')
args = parser.parse_args()

#######################################################
## extract biomart info to dict and store it to json ##
#######################################################
go_db = '/home/public/database/go/gene_ontology.1_2.obo'
go_anno_dict = {}
if args.biomart:
    if args.biomart.endswith('json'):
        with open(args.biomart) as biomart_info:
            go_anno_dict = json.load(biomart_info)
    else:
        reader = csv.reader(file(args.biomart, 'rb'))
        for n,each_record in enumerate(reader):
            if n!= 0:
                gene_id = each_record[0]
                go_id = each_record[1]
                go_domain = each_record[2]                
                go_des = each_record[3]
                gene_inter_id = each_record[4]
                gene_inter_des = each_record[5]
                if gene_id not in go_anno_dict:
                    go_anno_dict[gene_id] = [[go_id],[go_domain],[go_des],[gene_inter_id],[gene_inter_des]]
                else:
                    if go_id not in go_anno_dict[gene_id][0]:
                        go_anno_dict[gene_id][0].append(go_id)
                        go_anno_dict[gene_id][1].append(go_domain)                        
                        go_anno_dict[gene_id][2].append(go_des)
                    if gene_inter_id not in go_anno_dict[gene_id][2]:
                        go_anno_dict[gene_id][3].append(gene_inter_id)
                        go_anno_dict[gene_id][4].append(gene_inter_des)
        json_file = '%s.json' % args.biomart
        with open(json_file,'w') as json_file_info:
            json.dump(go_anno_dict, json_file_info)
elif args.go:
    if args.go.endswith('json'):
        with open(args.go) as go_info:
            go_anno_dict = json.load(go_info)
    else:
        go_db_info = obo_parser.GODag(go_db)        
        reader = csv.reader(file(args.go, 'rb'))
        for n,each_record in enumerate(reader):
            if n!= 0:
                gene_id = each_record[0]
                go_id = each_record[1]
                go_domain = go_db_info[go_id].namespace
                go_des = go_db_info[go_id].name
                if gene_id not in go_anno_dict:
                    go_anno_dict[gene_id] = [[go_id],[go_domain],[go_des]]
                else:
                    if go_id not in go_anno_dict[gene_id][0]:
                        go_anno_dict[gene_id][0].append(go_id)
                        go_anno_dict[gene_id][1].append(go_domain)                        
                        go_anno_dict[gene_id][2].append(go_des)                    
        go_json_file = '%s.json' % args.go
        with open(go_json_file,'w') as go_json_file_info:
            json.dump(go_anno_dict, go_json_file_info)
else:
    sys.exit('one of argument --biomart or --go is required')

############################
## extract gtf infomation ##
############################
tr_info_dict = {}
gene_name_dict = {}
with open(args.gtf) as gtf_info:
    if args.gtf.endswith('json'):
        tr_info_dict = json.load(gtf_info)
    else:
        tr_info_dict = RNAseq_tools.get_transcript_info(args.gtf)
        gtf_json_file = '%s.json' % args.gtf
        with open(gtf_json_file,'w') as gtf_json_file_info:
            json.dump(tr_info_dict,gtf_json_file_info)
for each_tr in tr_info_dict:
    if tr_info_dict[each_tr]['gene_name'] != '--':
        gene_id = tr_info_dict[each_tr]['gene_id']
        if gene_id not in gene_name_dict:
            gene_name_dict[gene_id] = tr_info_dict[each_tr]['gene_name']

###################################################
## extract KOBAS annotation and store it to json ##
###################################################
ko_anno_dict = {}
with open(args.kobas) as kobas_info:
    if args.kobas.endswith('json'):
        ko_anno_dict = json.load(kobas_info)
    else:
        tr_id = ''
        gene_id = ''
        for eachline in kobas_info:
            eachline_info = eachline.strip().split('\t')
            if len(eachline_info) == 2 and eachline_info[0].strip() == 'Query:':
                query_id = eachline_info[1]
                if args.type == 'gene':
                    query_id = tr_info_dict[query_id]['gene_id']
            if len(eachline_info) == 4 and eachline_info[2] == 'KEGG PATHWAY':
                ko_id = eachline_info[3]
                ko_des = eachline_info[1]
                if query_id not in ko_anno_dict:
                    ko_anno_dict[query_id] = [[ko_id], [ko_des]]
                else:
                    if ko_id not in ko_anno_dict[query_id][0]:
                        ko_anno_dict[query_id][0].append(ko_id)
                        ko_anno_dict[query_id][1].append(ko_des)
        ko_json_file = '%s.json' % args.kobas
        with open(ko_json_file,'w') as ko_json_file_info:
            json.dump(ko_anno_dict, ko_json_file_info)                    

#########################
## annotate diff table ##
#########################
output_info_list = []
with open(args.diff_out) as diff_out_info:
    for n, eachline in enumerate(diff_out_info):
        eachline_info = eachline.strip().split('\t') 
        gene_id = eachline_info[0] 
        if gene_name_dict:
            if n == 0:
                gene_name = 'Gene_name'
            else:
                gene_name = gene_name_dict[eachline_info[0]]
            eachline_info.insert(1,gene_name)        
                          
        if n == 0:
            if args.go:
                go_info = ['GO_term_accession', 'GO_domain', 'GO_term_name']
            else:                
                go_info = ['GO_term_accession', 'GO_domain', 'GO_term_name', 'InterPro_ID', 'InterPro_description']
        else:
            if args.go:
                go_info = ['--']*3
            else:
                go_info = ['--']*5            
            if gene_id in go_anno_dict:
                go_info = []
                for each in go_anno_dict[gene_id]:
                    go_info.append(','.join(each))
        eachline_info.extend(go_info)

        if n == 0:
            ko_info = ['KO_ID', 'KO_name']
        else:
            ko_info = ['--']*2
            if gene_id in ko_anno_dict:
                ko_info = []
                for each in ko_anno_dict[gene_id]:
                    ko_info.append(','.join(each))
        eachline_info.extend(ko_info)
        output_info_list.append('\t'.join(eachline_info))
python_tools.write_obj_to_file(output_info_list,args.output)










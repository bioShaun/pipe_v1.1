#! /usr/bin/python

'''
data:2016-3-9

1. Extract differentially expressed gene id and transform ensembl id to ncbi id.
2. Perform GO and KEGG enrichment analysis.

'''

import sys
import os
import argparse
import glob
LIB_LNCRNA = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, LIB_LNCRNA)
import circ
import json
import gzip
from HTSeq import GFF_Reader

mycwd = os.getcwd()

GENEID_DATA = '/home/public/gene_info'
GO_DATA = '/home/public/gene2go.gz'
GENEID_DATA_JSON = '/home/public/id_map.json'
GO_DATA_JSON = '/home/public/go_map.json'

RUN_GOSEQ = 'Rscript /home/lxgui/scripts/goseq.R'

parser = argparse.ArgumentParser()
parser.add_argument('--diff_out',help = 'DESeq2 differential analysis out directory.', required = True)
parser.add_argument('--compare', help = 'Compare name list.', required =True)
parser.add_argument('--qvalue', type = float ,help = 'qvalue cutoff.', default = 0.001)
parser.add_argument('--fc', type = float , help = 'log2foldchange cutoff.',default = 2.0)
parser.add_argument('--species' , help = 'KEGG species.', required =True)
parser.add_argument('--gtf' , help = 'GTF file.', required =True)
parser.add_argument('--go' , help = 'GO file, seperated with (,).', default = '')
parser.add_argument('--gene_length' , help = 'gene_length file, seperated with tab.', default = '')
parser.add_argument('--notransfer', action='store_true' ,help = 'not transfer ensembl id to NCBI id.', default = False)
parser.add_argument('--out_dir', help = 'output directory.', default = mycwd)
args = parser.parse_args()

compare_list = [each.strip() for each in open(args.compare)]
ensembl_to_ncbi_dict = {}
id_go_dict = {}

def diff_version(diff_line) :
    if len(diff_line) == 9 :
        return 'old'
    elif len(diff_line) == 7 :
        return 'new'
    else :
        return 'not DESeq2 result'

def get_go_map_dict() :
    global id_go_dict
    if not os.path.isfile(GO_DATA_JSON) :
        if not os.path.isfile(GO_DATA) :
            sys.exit('%s not exist!' % GO_DATA)
        with gzip.open(GO_DATA, "rb") as go_info :
            for n,eachline in enumerate(go_info) :
                if n != 0 :
                    each_go = eachline.strip().split('\t')
                    gene_id = each_go[1]
                    gene_go = each_go[2]
                    if gene_id not in id_go_dict :
                        id_go_dict[gene_id] = [gene_go]
                    else :
                        id_go_dict[gene_id].append(gene_go)
        with open(GO_DATA_JSON,'w') as go_map_file :                       
            json.dump(id_go_dict,go_map_file) 
    else :
        with open(GO_DATA_JSON,'r') as go_map_file :
            id_go_dict = json.load(go_map_file)                    

def get_id_map_dict():
    global ensembl_to_ncbi_dict
    if not os.path.isfile(GENEID_DATA_JSON) :
        if not os.path.isfile(GENEID_DATA) :
            sys.exit('%s not exist!' % GENEID_DATA)
        with open(GENEID_DATA,'r') as gene_info :
            for n,eachline in enumerate(gene_info) :
                if n != 0 :
                    each_info = eachline.strip().split('\t')
                    ncbi_id = each_info[1]                    
                    ensembl_id = each_info[2].upper()
                    if ensembl_id :
                        ensembl_to_ncbi_dict[ensembl_id] = ncbi_id
        with open(GENEID_DATA_JSON,'w') as id_map_file :                       
            json.dump(ensembl_to_ncbi_dict,id_map_file)            
    else :
        with open(GENEID_DATA_JSON,'r') as id_map_file :
            ensembl_to_ncbi_dict = json.load(id_map_file)

def get_diff_gene(diff_line,version_msg):
    gene_id = diff_line[0]
    if version_msg == 'old' :
        gene_fc = diff_line[4]
        gene_fdr = diff_line[8]
    elif version_msg == 'new' :
        gene_fc = diff_line[2]
        gene_fdr = diff_line[6]
    else :
        sys.exit(version_msg)  
    if gene_fc == 'NA' or gene_fdr == 'NA' :
        return None,None
    elif abs(float(gene_fc)) >= args.fc and float(gene_fdr) <= args.qvalue :
        ncbi_id = ''
        if gene_id in ensembl_to_ncbi_dict :
            ncbi_id = ensembl_to_ncbi_dict[gene_id]
        return gene_id,ncbi_id      
    else :
        return None,None

def get_genelist(diff_file,compare):
    diff_gene_list = []
    ncbi_gene_list = []
    for n,eachline in enumerate(open(diff_file,'r')) :            
        each_gene_exp = eachline.strip().split('\t')
        version_msg = diff_version(each_gene_exp)
        if n != 0 :
            diff_gene,ncbi_id = get_diff_gene(each_gene_exp,version_msg)
            if diff_gene :
                diff_gene_list.append(diff_gene)
            if ncbi_id :
                ncbi_gene_list.append(ncbi_id)
    return '\n'.join(diff_gene_list) + '\n','\n'.join(ncbi_gene_list) + '\n'

def get_transcript_length(gtf):
    tr_length_dict = {}
    for eachline in GFF_Reader(gtf):
        gene_id = eachline.attr['gene_id']
        transcript_id  = eachline.attr['transcript_id']
        start  = eachline.iv.start
        end    = eachline.iv.end
        length = end - start + 1
        if gene_id in tr_length_dict and transcript_id in tr_length_dict[gene_id] :
            tr_length_dict[gene_id][transcript_id] += length
        else :
            tr_length_dict.setdefault(gene_id,{})[transcript_id] = length
    return tr_length_dict

def Median(a):
    a = sorted(a)
    l = len(a)
    l2 = l / 2
    return a[l2] if l % 2 else (a[l2] + a[l2 - 1]) / 2.0

def get_gene_len(tr_length_dict,species,out_dir,notransfer = False) :
    gene_length_file = os.path.join(out_dir,'%s.gene_length' % species)
    with open(gene_length_file,'w') as gene_length_info :
        for each_gene in tr_length_dict :
            gene_length_list = []
            for each_tr in tr_length_dict[each_gene] :
                gene_length_list.append(tr_length_dict[each_gene][each_tr])
            each_gene_length = Median(gene_length_list)
            if notransfer :
                output_id = each_gene
                gene_length_info.write('%s\t%s\n' % (output_id,each_gene_length))        
            else :
                if each_gene in ensembl_to_ncbi_dict :
                    output_id = ensembl_to_ncbi_dict[each_gene]
                    gene_length_info.write('%s\t%s\n' % (output_id,each_gene_length))  
    return gene_length_file      
                

def get_go(tr_length_dict,species,out_dir,notransfer = False) :
    go_file = os.path.join(out_dir,'%s.go' % species)
    with open(go_file,'w') as go_info :
        go_info.write('GeneID,GO_accession\n')
        for each_gene in tr_length_dict :
            if notransfer :
                ncbi_id = each_gene
            else :
                if each_gene in ensembl_to_ncbi_dict :
                    ncbi_id = ensembl_to_ncbi_dict[each_gene]
                else :
                    continue
            if ncbi_id in id_go_dict :
                for each_go in id_go_dict[ncbi_id] :
                    go_info.write('%s,%s\n' % (ncbi_id,each_go))
    return go_file

def KEGG_enrich(infile,sp,output,format='id:ncbigene'):
    os.system('run_kobas.py -i {infile}  -t {format} -s {sp} -d K -o {output}'.format(**locals()))

def GO_enrich(difflist,length_file,go_file,output) :
    rungo_seq = RUN_GOSEQ
    os.system('{rungo_seq} {difflist} {length_file} {go_file} {output}'.format(**locals()))

def main() :

    # get_id_map_dict()
    # get_go_map_dict()
    # tr_length_dict = get_transcript_length(args.gtf)
    ## Check gene length and go file
    gene_length_file = args.gene_length
    if not args.gene_length :
        gene_length_file = get_gene_len(tr_length_dict,args.species,args.out_dir,args.notransfer)  
    go_file = args.go
    if not go_file :
        go_file = get_go(tr_length_dict,args.species,args.out_dir,args.notransfer)

    for each_compare in compare_list :
        # ## DESeq2 result
        # each_compare_file = glob.glob(r'%s/*%s.DESeq2.DE_results' % (args.diff_out,each_compare))[0]
        # ## Extract Diff id list and transform to NCBI id
        # diff_genes,ncib_genes = get_genelist(each_compare_file,each_compare)        
        # out_file1 = os.path.join(args.out_dir,'%s.diffgene.list' % each_compare)    
        out_file2 = os.path.join(args.out_dir,'ncbi.%s.diffgene.list' % each_compare)    
        # circ.write_str_to_file(diff_genes,out_file1)
        # circ.write_str_to_file(ncib_genes,out_file2)
        # ## KEGG enrichment 
        # kegg_out = os.path.join(args.out_dir,'%s.KEGG.enrich.xls' % each_compare) 
        # KEGG_enrich(out_file2,args.species,kegg_out)  
        ## GO enrichment
        go_out = os.path.join(args.out_dir,'%s.GO.enrich.xls' % each_compare)
        GO_enrich(out_file2, gene_length_file, go_file, go_out)
        
if __name__ == '__main__':
    main()
                



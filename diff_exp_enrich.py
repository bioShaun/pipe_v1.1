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
KO_PEP_DIR = '/home/lxgui/test_kobas/seq_pep/'

RUN_GOSEQ = 'Rscript /home/lxgui/scripts/goseq.R'
BLAST_BIN = '/usr/bin/'

parser = argparse.ArgumentParser()
parser.add_argument('--diff_out',help = 'DESeq2 differential analysis out directory.', required = True)
parser.add_argument('--compare', help = 'Compare name list.', required =True)
parser.add_argument('--qvalue', type = float ,help = 'qvalue cutoff.', default = 0.001)
parser.add_argument('--fc', type = float , help = 'log2foldchange cutoff.',default = 2.0)
parser.add_argument('--species' , help = 'KEGG species.', required =True)
parser.add_argument('--gtf' , help = 'GTF file.', required =True)
parser.add_argument('--go' , help = 'GO file, seperated with (,).', default = '')
parser.add_argument('--seq' , help = '[cds.fa,pep.fa]', default = '')
parser.add_argument('--blast_out' , help = 'blast out file', default = '')
parser.add_argument('--gene_length' , help = 'gene_length file, seperated with tab.', default = '')
parser.add_argument('--notransfer', action='store_true' ,help = 'not transfer ensembl id to NCBI id.', default = False)
parser.add_argument('--nokegg', action='store_true' ,help = 'not perform kegg analysis.', default = False)
parser.add_argument('--nogo', action='store_true' ,help = 'not extract go analysis.', default = False)
parser.add_argument('--out_dir', help = 'output directory.', default = mycwd)
args = parser.parse_args()

compare_list = [each.strip() for each in open(args.compare)]
ensembl_to_ncbi_dict = {}
id_go_dict = {}

# def get_stat_so():
#     jobs_status = os.popen('ps -ux').read()
#     if 'javareconf' not in jobs_status :
#         os.system('sh ~/.R_java.sh')

def diff_version(diff_line) :
    if len(diff_line) == 9 :
        return 'old_DESeq2'
    elif len(diff_line) == 6 :
        return 'new_DESeq2'
    elif len(diff_line) == 4 :
        return 'edgeR'
    elif len(diff_line) == 8 :
        return 'Cuffdiff'    
    else :
        return 0

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
    if version_msg == 'old_DESeq2' :
        gene_fc = diff_line[4]
        gene_fdr = diff_line[8]
    elif version_msg == 'new_DESeq2' :
        gene_fc = diff_line[2]
        gene_fdr = diff_line[6]
    elif version_msg == 'edgeR' :
        gene_fc = diff_line[1]
        gene_fdr = diff_line[4]   
    elif version_msg == 'Cuffdiff' :
        gene_fc = diff_line[5]
        gene_fdr = diff_line[7]         
    else :
        sys.exit('not DEseq or edgeR result')  
    if gene_fc == 'NA' or gene_fdr == 'NA' :
        return None
    elif abs(float(gene_fc)) >= args.fc and float(gene_fdr) <= args.qvalue :
        return gene_id
    else :
        return None

def get_genelist(diff_file,compare):
    diff_gene_list = []
    for n,eachline in enumerate(open(diff_file,'r')) :            
        each_gene_exp = eachline.strip().split('\t')
        if n == 0 :
            version_msg = diff_version(each_gene_exp)
        else :
            diff_gene = get_diff_gene(each_gene_exp,version_msg)
            if diff_gene :
                diff_gene_list.append(diff_gene)
    return diff_gene_list

def id_list_convertion(ensembl_id) :
    ncbi_id_list = []
    for each in ensembl_id :
        if each in ensembl_to_ncbi_dict :
            ncbi_id_list.append(ensembl_to_ncbi_dict[each])
    return ncbi_id_list

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

def check_blast_program(seq_file):
    seq_file_name = os.path.basename(seq_file)
    if 'pep' in seq_file_name :
        blast_program = 'blastp'
    elif 'cds' in seq_file_name :
        blast_program = 'blastx'
    else :
        sys.exit('Do not know use which blast program !')
    return blast_program

def check_blast_database(seq_file) :
    blast_pr = check_blast_program(seq_file)
    suffix_list = ['phr','pin','psq']
    for each in suffix_list :
        if not os.path.isfile('%s.%s' % (seq_file,each)) :
            os.system('%s/formatdb -i %s -p T' % (BLAST_BIN,seq_file))
            break

def get_blast_out(seq_file,species,out_dir):
    blast_bin = BLAST_BIN
    ko_seq = os.path.join(KO_PEP_DIR,'%s.pep.fasta' % species)
    check_blast_database(ko_seq)
    blast_pr = check_blast_program(seq_file)
    os.system('{blast_bin}/{blast_pr}  -query {seq_file} -db {ko_seq} -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -num_threads 1 -out {out_dir}/{species}.blasttab'.format(**locals()))
    print '{blast_bin}/{blast_pr}  -query {seq_file} -db {ko_seq} -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -num_threads 1 -out {out_dir}/{species}.blasttab'.format(**locals())
    return '{out_dir}/{species}.blasttab'.format(**locals())

def convert_blast(gtf,kegg_blast,out_dir,species):
    tr_gene_dict = {}
    for eachline in GFF_Reader(gtf):
        gene_id = eachline.attr['gene_id']
        transcript_id  = eachline.attr['transcript_id']     
        tr_gene_dict[transcript_id] = gene_id
    diff_blast_dict = {}
    with open(kegg_blast,'r') as kegg_blast_info :
        for eachline in kegg_blast_info :
            each_tr_blast = eachline.strip().split('\t')
            gene_id = tr_gene_dict[each_tr_blast[0]]
            if gene_id not in diff_blast_dict :
                diff_blast_dict[gene_id] = each_tr_blast[1:]
            else :
                if diff_blast_dict[gene_id][-2] == each_tr_blast[-2] and diff_blast_dict[gene_id][-1] < each_tr_blast[-1]:
                    diff_blast_dict[gene_id] = each_tr_blast[1:]
                elif diff_blast_dict[gene_id][-2] > each_tr_blast[-2] :
                    diff_blast_dict[gene_id] = each_tr_blast[1:]
                else :
                    continue
    diff_blast = os.path.join(out_dir,'%s.gene.blasttab.json' % species)
    with open(diff_blast,'w') as diff_blast_file :
        json.dump(diff_blast_dict,diff_blast_file)
    return diff_blast

def extract_diff_blast(diff_list,gene_blast_json,compare,out_dir):
    with open(gene_blast_json,'r') as gene_blast_info :
        gene_blast_dict = json.load(gene_blast_info)
    diff_gene_list = [each.strip() for each in open(diff_list,'r')]
    diff_blast_out_file = os.path.join(out_dir,'%s.blasttab' % compare)
    with open(diff_blast_out_file,'w') as diff_blast_out :
        for each_gene in diff_gene_list :
            if each_gene in gene_blast_dict :
                blast_info = '\t'.join(gene_blast_dict[each_gene])
                diff_blast_out.write('%s\t%s\n' % (each_gene,blast_info))
    return diff_blast_out_file

def KEGG_enrich(infile,sp,output):
    if infile.endswith('.list') :
        os.system('run_kobas.py -i {infile}  -t id:ncbigene -s {sp} -d K -o {output}'.format(**locals()))
    elif infile.endswith('.blasttab') :
        os.system('run_kobas.py -i {infile}  -t blastout:tab -s {sp} -d K -o {output}'.format(**locals()))

def GO_enrich(difflist,length_file,go_file,output) :
    rungo_seq = RUN_GOSEQ
    os.system('{rungo_seq} {difflist} {length_file} {go_file} {output}'.format(**locals()))

def summery_enrich_results(out_dir):
    working_directory = os.path.join(args.out_dir,'processing')
    diff_list_dir = os.path.join(out_dir,'diffgene_list')
    go_results = os.path.join(out_dir,'go_enrichment')
    kegg_results = os.path.join(out_dir,'kegg_enrichment')
    map(circ.circ_mkdir_unix,[diff_list_dir,go_results,kegg_results])
    os.system('cp %s/*diffgene.list %s' % (working_directory,diff_list_dir))
    os.system('cp %s/*GO.enrich.xls %s' % (working_directory,go_results))
    os.system('cp %s/*KEGG.enrich.xls %s' % (working_directory,kegg_results))

def main() :
    working_directory = os.path.join(args.out_dir,'processing')
    if not os.path.exists(working_directory) :
        circ.circ_mkdir_unix(working_directory)
    # get_stat_so()
    gene_length_file = args.gene_length
    go_file = args.go
    tr_length_dict = {}
    ## Prepare gene length and go file
    if not args.gene_length :
        tr_length_dict = get_transcript_length(args.gtf)
        gene_length_file = get_gene_len(tr_length_dict,args.species,working_directory,args.notransfer)  
    if not args.go :
        if not tr_length_dict :
            tr_length_dict = get_transcript_length(args.gtf)
        go_file = get_go(tr_length_dict,args.species,working_directory,args.notransfer)        

    ko_blast = args.blast_out
    if args.seq :
        ## KEGG blast
        ko_blast = get_blast_out(args.seq,args.species,working_directory)

    for each_compare in compare_list :
        ## Diff result
        each_compare_file = glob.glob(r'%s/*%s.*.DE_results' % (args.diff_out,each_compare))[0]
        ## Extract Diff id list and transform to NCBI id        
        diff_list_file = os.path.join(working_directory,'%s.diffgene.list' % each_compare)
        diff_genes = get_genelist(each_compare_file,each_compare)          
        circ.write_obj_to_file(diff_genes,diff_list_file)
        if not args.notransfer :
            diff_list_file = os.path.join(working_directory,'ncbi.%s.diffgene.list' % each_compare)
            diff_genes = id_list_convertion(diff_genes)
            circ.write_obj_to_file(diff_genes,diff_list_file)

        ## KEGG enrichment
        if not args.nokegg :
            kegg_out = os.path.join(working_directory,'%s.KEGG.enrich.xls' % each_compare) 
            if not args.notransfer :
                KEGG_enrich(diff_list_file,args.species,kegg_out)  
            else :
                if '.json' not in ko_blast :
                    ko_blast = convert_blast(args.gtf,ko_blast,working_directory,args.species)
                diff_blast_out = extract_diff_blast(diff_list_file,ko_blast,each_compare,working_directory)
                KEGG_enrich(diff_blast_out,args.species,kegg_out)

        # GO enrichment
        if not args.nogo :
            go_out = os.path.join(working_directory,'%s.GO.enrich.xls' % each_compare)
            GO_enrich(diff_list_file, gene_length_file, go_file, go_out)

    ## enrichment analysis results summary
    summery_enrich_results(args.out_dir)

if __name__ == '__main__':
    main()
                



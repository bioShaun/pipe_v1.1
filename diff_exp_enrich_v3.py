'''
data:2016-3-24

1. Extract differentially expressed gene id and transform ensembl id to ncbi id.
2. Perform GO and KEGG enrichment analysis.
3. Add up and down gene list enrichment analysis
4. Add enrichment plot

last updat:2016-03-30
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
import RNAseq_tools

mycwd = os.getcwd()

GENEID_DATA = '/home/public/gene_info'
GO_DATA = '/home/public/gene2go.gz'
GENEID_DATA_JSON = '/home/public/id_map.json'
GO_DATA_JSON = '/home/public/go_map.json'
KO_PEP_DIR = '/home/lxgui/test_kobas/seq_pep/'

RUN_GOSEQ = 'Rscript /home/lxgui/scripts/goseq.R'
BLAST_BIN = '/usr/bin/'
PLOT_R_BIN = '/home/public/R_packages/R-3.1.0/bin/'
SIMPLE_PLOT_SP = '/home/lxgui/scripts/enrich_bar.2016-3-15.R'

parser = argparse.ArgumentParser()
parser.add_argument('--diff_out',help = 'DESeq2 differential analysis out directory.', required = True)
parser.add_argument('--diff_type',help = 'gene or transcript differential analysis out.',choices = ['gene','transcript'], default = 'gene')
parser.add_argument('--compare', help = 'Compare name list.', required =True)
parser.add_argument('--qvalue', type = float ,help = 'qvalue cutoff.', default = 0.001)
parser.add_argument('--fc', type = float , help = 'log2foldchange cutoff.',default = 2.0)
parser.add_argument('--species' , help = 'KEGG species.', required =True)
parser.add_argument('--gtf' , help = 'GTF file.', required =True)
parser.add_argument('--go' , help = 'GO file, seperated with (,).', default = '')
parser.add_argument('--seq' , help = '[cds.fa,pep.fa]', default = '')
parser.add_argument('--blast_out' , help = 'blast out file', default = '')
parser.add_argument('--gene_length' , help = 'gene_length file, seperated with tab.', default = '')
parser.add_argument('--nokegg', action='store_true' ,help = 'not perform kegg analysis.', default = False)
parser.add_argument('--nogo', action='store_true' ,help = 'not extract go analysis.', default = False)
parser.add_argument('--out_dir', help = 'output directory.', default = mycwd)
args = parser.parse_args()

compare_list = [each.strip() for each in open(args.compare)]
id_go_dict = {}

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
        sys.exit('not DEseq or edgeR or Cuffdiff result')  
    if gene_fc == 'NA' or gene_fdr == 'NA' :
        return None,None
    elif abs(float(gene_fc)) >= args.fc and float(gene_fdr) <= args.qvalue :
        if float(gene_fc) > 0 :
            return gene_id,'up'
        else :
            return gene_id,'down'
    else :
        return None,None

def get_genelist(diff_file,diff_gene_dict):
    for n,eachline in enumerate(open(diff_file,'r')) :            
        each_gene_exp = eachline.strip().split('\t')
        if n == 0 :
            version_msg = diff_version(each_gene_exp)
        else :
            diff_gene,reg_flag = get_diff_gene(each_gene_exp,version_msg)
            if diff_gene :
                if 'all' in diff_gene_dict :
                    diff_gene_dict['all'].append(diff_gene) 
                else :
                    diff_gene_dict['all'] = [diff_gene]
                if reg_flag in diff_gene_dict :
                    diff_gene_dict[reg_flag].append(diff_gene)
                else :
                    diff_gene_dict[reg_flag] = [diff_gene]

def summery_enrich_results(out_dir):
    working_directory = os.path.join(args.out_dir,'processing')
    diff_list_dir = os.path.join(out_dir,'diffgene_list')
    go_results = os.path.join(out_dir,'go_enrichment')
    kegg_results = os.path.join(out_dir,'kegg_enrichment')
    map(circ.circ_mkdir_unix,[diff_list_dir,go_results,kegg_results])
    os.system('cp %s/*diffgene.list %s' % (working_directory,diff_list_dir))
    os.system('cp %s/*GO.enrich.xls %s' % (working_directory,go_results))
    os.system('cp %s/*GO.bar.* %s' % (working_directory,go_results))
    os.system('cp %s/*KEGG.enrich.xls %s' % (working_directory,kegg_results))
    os.system('cp %s/*KEGG.bar.* %s' % (working_directory,kegg_results))

if __name__ == '__main__':
    reg_list = ['all','up','down']
    working_directory = os.path.join(args.out_dir,'processing')
    if not os.path.exists(working_directory) :
        circ.circ_mkdir_unix(working_directory)
    # get_stat_so()
    
    ## Prepare gene length and go file     
    my_go_enrich = RNAseq_tools.GO_enrich()
    my_go_enrich.target_length = args.gene_length
    my_go_enrich.gtf = args.gtf  
    my_go_enrich.go = args.go
    my_go_enrich.out_dir = working_directory      
    my_go_enrich.target_type = args.diff_type
    my_go_enrich.get_target_length_table()
    if args.nogo :
        my_go_enrich.auot_run = False

    ## KEGG annotate
    my_kegg_enrich = RNAseq_tools.KEGG_enrich()
    my_kegg_enrich.target_type = args.diff_type 
    my_kegg_enrich.species = args.species          
    my_kegg_enrich.out_dir = working_directory
    if args.nokegg :
        my_kegg_enrich.auot_run = False
    if not args.gene_length :
        my_kegg_enrich.transcript_dict = my_go_enrich.transcript_dict
    else :
        my_kegg_enrich.gtf = args.gtf
    if args.blast_out :
        my_kegg_enrich.all_blast_out = args.blast_out            
    elif args.seq :
        my_kegg_enrich.seq = args.seq                   
        my_kegg_enrich.all_blast_out = os.path.join(working_directory,'%s.all.blasttab' % args.species)
        my_kegg_enrich.get_blast_out()
    else :
        sys.exit('you need to input either species cds/pep blast results or sequence.')
    my_kegg_enrich.convert_blast() 

    for each_compare in compare_list :
        each_compare_outdir = os.path.join(working_directory,each_compare)
        circ.circ_mkdir_unix(each_compare_outdir)
        ## Diff result
        each_compare_file = glob.glob(r'%s/*%s.*.DE_results' % (args.diff_out,each_compare))[0]
        ## Extract Diff id list 
        each_diff_dict = {}
        get_genelist(each_compare_file,each_diff_dict)  
        for each_reg in reg_list : 
            if each_reg not in each_diff_dict :
                continue
            diff_list = each_diff_dict[each_reg]
            diff_list_file = os.path.join(each_compare_outdir,'%s.%s.diffgene.list' % (each_compare,each_reg))
            circ.write_obj_to_file(diff_list,diff_list_file)
            ## add KEGG enrichment info
            my_kegg_enrich.add_enrich_info(each_compare,each_reg,diff_list)     
            ## add GO enrichment info       
            my_go_enrich.add_enrich_info(each_compare,each_reg,diff_list_file)

    ## run kegg enrichment
    if not args.nokegg :
        my_kegg_enrich.run_KEGG_enrich()
    ## run go enrichment
    if not args.nogo :
        my_go_enrich.run_GO_enrich()

    # ## enrichment analysis results summary
    # summery_enrich_results(args.out_dir)
                



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
import python_tools

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
parser.add_argument('--species' , help = 'KEGG species.', required =True)
parser.add_argument('--background' , help = 'KEGG analysis backgroud.', default = "")
parser.add_argument('--diffmethod', help = 'differential analysis software', required = True)
parser.add_argument('--gtf' , help = 'GTF file.', required =True)
parser.add_argument('--go' , help = 'GO file, seperated with (,).', default = '')
parser.add_argument('--topgo_go', help = 'gene go map file', default = '')
parser.add_argument('--seq' , help = '[cds.fa,pep.fa]', default = '')
parser.add_argument('--blast_out' , help = 'blast out file', default = '')
parser.add_argument('--gene_length' , help = 'gene_length file, seperated with tab.', default = '')
parser.add_argument('--nokegg', action='store_true' ,help = 'not perform kegg analysis.', default = False)
parser.add_argument('--nogo', action='store_true' ,help = 'not extract go analysis.', default = False)
parser.add_argument('--kegg_pathview', action='store_true', help = 'run pathview analysis', default = False)
parser.add_argument('--topgo', action='store_true', help = 'run topgo DAG plot', default = False)
parser.add_argument('--sRNA', action='store_true', help = 'sRNA analysis', default = False)
parser.add_argument('--sRNA_target', help = 'sRNA target pairs', default = '')
parser.add_argument('--out_dir', help = 'output directory.', default = mycwd)
args = parser.parse_args()

if not args.background:
    kegg_background = args.species
else:
    kegg_background = args.background

compare_list = [each.strip() for each in open(args.compare)]
id_go_dict = {}

def get_diff_gene(difffile,outgene):
    diff_gene_list = []
    if os.path.isfile(difffile):
        with open(difffile) as difffile_info:
            for n,eachline in enumerate(difffile_info):
                if n!= 0:
                    eachline_info = eachline.split('\t')
                    diff_gene_list.append(eachline_info[0])
    else:
        sys.exit('%s not exists!' % difffile)
    python_tools.write_obj_to_file(diff_gene_list,outgene)

def get_diff_target(difffile, sRNA_target_dict, out_target):
    target_mRNA_list = []
    diff_list = [each.strip() for each in open(difffile)]
    for each_miRNA in diff_list:
        if each_miRNA in sRNA_target_dict:            
            target_mRNA_list.extend(sRNA_target_dict[each_miRNA])
    target_mRNA_list = list(set(target_mRNA_list))
    python_tools.write_obj_to_file(target_mRNA_list, out_target)

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

    ## KEGG annotate
    my_kegg_enrich = RNAseq_tools.KEGG_enrich()
    if not args.nokegg:
        my_kegg_enrich.target_type = args.diff_type 
        my_kegg_enrich.species = args.species          
        my_kegg_enrich.out_dir = working_directory
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

    ## prepare miRNA anaylis
    if args.sRNA:
        sRNA_target = args.sRNA_target
        sRNA_target_dict = {}
        with open(sRNA_target) as sRNA_target:
            for eachline in sRNA_target:
                eachline_info = eachline.strip().split('\t')
                miRNA = eachline_info[0]
                mRNA_gene = eachline_info[1]
                # mRNA_gene = my_kegg_enrich.transcript_dict[mRNA]['gene_id']
                if miRNA not in sRNA_target_dict:
                    sRNA_target_dict[miRNA] = [mRNA_gene]
                else:
                    sRNA_target_dict[miRNA].append(mRNA_gene)

    for each_compare in compare_list :
        each_compare_outdir = os.path.join(working_directory,each_compare)
        circ.circ_mkdir_unix(each_compare_outdir)
        ## Diff result
        cond_list = each_compare.split('_vs_')
        diff_gene_file_list = []
        for each_cond in cond_list:
            each_cond_up_file = glob.glob(r'%s/*%s.*.DE_results.*.%s-UP.subset' % (args.diff_out,each_compare,each_cond))[0]
            each_cond_up_gene = os.path.join(each_compare_outdir,'%s-UP.list' % each_cond)
            get_diff_gene(each_cond_up_file,each_cond_up_gene)
            if args.sRNA:
                each_cond_up_target = os.path.join(each_compare_outdir,'%s-UP-target.list' % each_cond)
                get_diff_target(each_cond_up_gene, sRNA_target_dict, each_cond_up_target)
                diff_gene_file_list.append(each_cond_up_target)
            else:
                diff_gene_file_list.append(each_cond_up_gene)
        if args.sRNA:                
            all_diff_gene = os.path.join(each_compare_outdir,'ALL-target.list')
        else:
            all_diff_gene = os.path.join(each_compare_outdir,'ALL.list')
        python_tools.merge_files(diff_gene_file_list,all_diff_gene)
        diff_gene_file_list.append(all_diff_gene)
        ## Extract Diff id list 
        for each_gene_list in diff_gene_file_list :
            each_name = os.path.basename(each_gene_list)
            each_name = os.path.splitext(each_name)[0]
            list_obj = [each.strip() for each in open(each_gene_list)]
            if list_obj:
                ## add KEGG enrichment info
                my_kegg_enrich.add_enrich_info(each_compare,each_name,list_obj)     
                ## add GO enrichment info       
                my_go_enrich.add_enrich_info(each_compare,each_name,each_gene_list)

    ## run kegg enrichment
    if not args.nokegg :
        my_kegg_enrich.run_KEGG_enrich_old(kegg_background)

    ## run go enrichment
    if not args.nogo :
        my_go_enrich.run_GO_enrich()

    if args.kegg_pathview:
        if args.sRNA:
            my_kegg_enrich.sRNA_target = args.sRNA_target
            my_kegg_enrich.run_kegg_pathview(args.diff_out, args.diffmethod, True)
        else:
            my_kegg_enrich.run_kegg_pathview(args.diff_out, args.diffmethod)

    if args.topgo:
        my_go_enrich.run_GO_DAG(args.topgo_go)

    # ## enrichment analysis results summary
    # summery_enrich_results(args.out_dir)
                



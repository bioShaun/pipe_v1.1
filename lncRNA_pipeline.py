'''
data:2016-06-03

lncRNA seq pipeline
1. mapping

'''

import sys
import os
import argparse
import RNAseq_tools
import python_tools

mycwd = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument('--proj_name', help = 'project name', required = True)
parser.add_argument('--fq_dir',help = 'clean fq file directory.', required = True)
parser.add_argument('--sample_list',help = 'samples for analysis.', default = "")
parser.add_argument('--genome_bowtie2_index',help = 'geneome bowtie2 index.', default = '')
parser.add_argument('--trans_bowtie2_index',help = 'transcriptome bowtie2 index.', default = '')
parser.add_argument('--gtf',help = 'gtf file path.', default = '')
parser.add_argument('--libtype',help = 'library type.', choices = ['fr-unstranded','fr-firststrand','fr-secondstrand'], default = 'fr-firststrand')
parser.add_argument('--mapping_dir',help = 'clean fq file directory.', default = '')
parser.add_argument('--mapping_thread',help = 'mapping threads.', default = 8)
parser.add_argument('--assembly_dir',help = 'assembly directory.', default = '')
parser.add_argument('--assembly_thread',help = 'assembly threads.', default = 8)
parser.add_argument('--group_sample',help = 'group sample mapping file.', default = "")
parser.add_argument('--gene_transcript', help = 'gene transcript id mapping file .', default = "")
parser.add_argument('--transcript_fa' , help = 'transcript fasta file.', default = "")
parser.add_argument('--gene_length' , help = 'gene length file.', default = "")
parser.add_argument('--go_anno' , help = 'gene go annotation file.', default = "")
parser.add_argument('--topgo_anno' , help = 'gene topgo annotation file.', default = "")
parser.add_argument('--kegg_blast' , help = 'kegg blast annotation file.', default = "")
parser.add_argument('--kegg_species' , help = 'kegg species abbr.', default = "")
parser.add_argument('--kegg_background' , help = 'kegg background species.', default = "")
parser.add_argument('--quant_method', help = 'quantification software', default = 'salmon')
parser.add_argument('--quant_dir', help = 'quantification output directory.', default = "")
parser.add_argument('--diff_dir', help = 'differential analysis output directory.', default = "")
parser.add_argument('--enrich_dir', help = 'enrichment analysis output directory.', default = "")
parser.add_argument('--platform', help = 'analysis platform.', default = 'th')
parser.add_argument('--mapping', help = 'run mapping', action='store_true', default = False)
parser.add_argument('--assembly', help = 'run assembly', action='store_true', default = False)
parser.add_argument('--quant', help = 'run quantification', action='store_true', default = False)
parser.add_argument('--go', help = 'run go enrichment analysis', action='store_true', default = False)
parser.add_argument('--kegg', help = 'run kegg enrichment analysis', action='store_true', default = False)
parser.add_argument('--quant_anno', help = 'quantification annotation files(sep with ",")', default = "")
parser.add_argument('--qc_dir', help = 'qc directory', default = '')
parser.add_argument('--qc', help = 'run qc', action='store_true', default = False)
parser.add_argument('--analysis', help = 'analysis pipeline ', choices = ['mRNA', 'lncRNA', 'miRNA'], default = 'mRNA')
args = parser.parse_args()

if __name__ == '__main__':
    my_project = RNAseq_tools.RNAseq_pipeline()
    ## parameter assignment
    my_project.cleandata_dir = args.fq_dir
    if args.sample_list:
        my_project.sample_list = [each.strip() for each in open(args.sample_list)]
    my_project.group_sample = args.group_sample
    my_project.gene_trans = args.gene_transcript
    my_project.transcript_fa = args.transcript_fa
    my_project.quant_program = args.quant_method
    my_project.quant_dir = args.quant_dir
    my_project.diff_dir = args.diff_dir
    my_project.platform = args.platform
    my_project.genome_bowtie2_index = args.genome_bowtie2_index
    my_project.trans_bowtie2_index = args.trans_bowtie2_index
    my_project.gtf = args.gtf
    my_project.libtype = args.libtype
    my_project.mapping_dir = args.mapping_dir
    my_project.mapping_thread = args.mapping_thread
    my_project.assembly_dir = args.assembly_dir
    my_project.assembly_thread = args.assembly_thread
    my_project.run_go = args.go
    my_project.run_kegg = args.kegg
    my_project.go = args.go_anno
    my_project.topgo = args.topgo_anno
    my_project.blast_out = args.kegg_blast
    my_project.gene_length = args.gene_length
    my_project.enrich_dir = args.enrich_dir
    my_project.kegg_species = args.kegg_species
    if args.kegg_background:
        my_project.kegg_background = args.kegg_background
    else:
        my_project.kegg_background = args.kegg_species
    my_project.anno_files = args.quant_anno
    my_project.qc_dir = args.qc_dir
    
#    my_project.qvalue = 0.001
#    my_project.logfc = 2
    ## qc
    if args.qc:
        my_project.run_qc()
    ## mapping
    if args.mapping:
        my_project.run_mapping()
    ## asssembly
    if args.assembly:
        my_project.run_assembly()
    ## quantification    
    if args.quant:
        my_project.run_quant()
    ## enrichment
    if args.go or args.kegg:
        my_project.run_enrich()
    ## result
    if args.analysis == 'mRNA':
        my_project.get_mRNA_results(mycwd, args.proj_name)

'''
get circos diff distribution plot data

'''

import sys
import os
import argparse
from HTSeq import GFF_Reader

parser = argparse.ArgumentParser()
parser.add_argument('--gtf', help = 'gtf file', required = True)
parser.add_argument('--diff_out', help = 'diff analysis result', required = True)
parser.add_argument('--out_dir', help = 'output directory', required = True)
parser.add_argument('--qvalue', type = float ,help = 'qvalue cutoff.', default = 0.001)
parser.add_argument('--fc', type = float , help = 'log2foldchange cutoff.',default = 2.0)
args = parser.parse_args()

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
        gene_fc = float(gene_fc)
        if float(gene_fc) > 0 :
            return gene_id,gene_fc
        else :
            return gene_id,gene_fc
    else :
        return None,None

def get_genelist(diff_file,diff_gene_dict):
    for n,eachline in enumerate(open(diff_file,'r')) :            
        each_gene_exp = eachline.strip().split('\t')
        if n == 0 :
            version_msg = diff_version(each_gene_exp)
        else :
            diff_gene,gene_fc = get_diff_gene(each_gene_exp,version_msg)
            if diff_gene :
                diff_gene_dict[diff_gene] = gene_fc

diff_dict = {}
get_genelist(args.diff_out,diff_dict)  
mRNA_diff_file = os.path.join(args.out_dir,'mRNA_diff_hist.txt')
ncRNA_diff_file = os.path.join(args.out_dir,'ncRNA_diff_hist.txt')
mRNA_diff_file_info = open(mRNA_diff_file,'w')
ncRNA_diff_file_info = open(ncRNA_diff_file,'w')

for eachline in GFF_Reader(args.gtf) :
    if eachline.type == 'gene' :
        gene_id = eachline.attr['gene_id']
        if gene_id in diff_dict :
            gene_type = eachline.attr['gene_biotype']
            start = eachline.iv.start
            end = eachline.iv.end
            chrom = 'chr%s' % eachline.iv.chrom
            gene_fc = diff_dict[gene_id]
            if gene_type == 'protein_coding' :
                mRNA_diff_file_info.write('%s %s %s %s\n' % (chrom,start,end,gene_fc))
            else :
                ncRNA_diff_file_info.write('%s %s %s %s\n' % (chrom,start,end,gene_fc))
mRNA_diff_file_info.close()
ncRNA_diff_file_info.close()



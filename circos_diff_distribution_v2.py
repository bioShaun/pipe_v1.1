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
        return None,None,None        
    elif abs(float(gene_fc)) >= args.fc and float(gene_fdr) <= args.qvalue :
        gene_fc = float(gene_fc)
        if float(gene_fc) > 0 :
            return gene_id,gene_fc,'up'
        else :
            return gene_id,gene_fc,'down'
    else :
        gene_fc = float(gene_fc)
        return gene_id,gene_fc,'no'


def get_genelist(diff_file,diff_gene_dict):
    for n,eachline in enumerate(open(diff_file,'r')) :            
        each_gene_exp = eachline.strip().split('\t')
        if n == 0 :
            version_msg = diff_version(each_gene_exp)
        else :
            diff_gene,gene_fc,reg = get_diff_gene(each_gene_exp,version_msg)
            if diff_gene :
                diff_gene_dict[diff_gene] = [gene_fc,reg]

def get_color(gene_type,fc,reg):
    if reg == 'no' :
        if gene_type == 'protein_coding' :
            if fc > 0 :
                return 'color=vvlred'
            else :
                return 'color=vvlblue'
        else :
            if fc > 0 :
                return 'color=vvlorange'
            else :
                return 'color=vvlgreen'            
    elif reg == 'up' :
        if gene_type == 'protein_coding' :
            return 'color=red'
        else :
            return 'color=orange'
    elif reg == 'down' :
        if gene_type == 'protein_coding' :
            return 'color=blue'
        else :
            return 'color=green'        
    else :
        sys.exit('wrong reg type : %s !' % reg)

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
            gene_fc = diff_dict[gene_id][0]
            reg = diff_dict[gene_id][1]
            color = get_color(gene_type,gene_fc,reg)
            if gene_type == 'protein_coding' :
                mRNA_diff_file_info.write('%s %s %s %s %s\n' % (chrom,start,end,gene_fc,color))
            else :
                ncRNA_diff_file_info.write('%s %s %s %s %s\n' % (chrom,start,end,gene_fc,color))
mRNA_diff_file_info.close()
ncRNA_diff_file_info.close()



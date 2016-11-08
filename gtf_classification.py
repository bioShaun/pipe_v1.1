import sys
import os
import argparse
from HTSeq import GFF_Reader
import json

parser = argparse.ArgumentParser(description = 'Filter transcripts by expression.')
parser.add_argument('--gtf', help = 'input gtf file to be classified.', required = True)
parser.add_argument('--out_dir', help = 'output file directory.' , required = True)
argv = vars(parser.parse_args())

gtf = argv['gtf'].strip()
out_dir = argv['out_dir'].strip()

mRNA_tag_list = 'mRNA,protein_coding'.split(',')
lncRNA_tag_list = 'lncRNA,3prime_overlapping_ncrna,ambiguous_orf,antisense,antisense_RNA,lincRNA,ncrna_host,processed_transcript,sense_intronic,sense_overlapping'.split(',')
lncRNA_candidate_tag_list = 'ncRNA'.split(',')
mRNA_gtf = os.path.join(out_dir,'mRNA.gtf')
lncRNA_gtf = os.path.join(out_dir,'lncRNA.gtf')
other_ncRNA_gtf = os.path.join(out_dir,'other_ncRNA.gtf')

def get_genetype_length(gffreader_line):
    line = gffreader_line.get_gff_line().strip()
    start = gffreader_line.iv.start
    end = gffreader_line.iv.end  
    exon_length = end - start + 1
    if 'gene_type' in gffreader_line.attr :
        gene_type = gffreader_line.attr['gene_type']
    elif 'gene_biotype' in gffreader_line.attr :
        gene_type = gffreader_line.attr['gene_biotype']
    else :
        sys.exit('no gene type in gtfline : {line}'.format(**locals()))
    return gene_type,exon_length

if __name__ == '__main__':
    candidate_lncRNA_dict = {}
    mRNA_gtf_file = open(mRNA_gtf,'w')
    lncRNA_gtf_file = open(lncRNA_gtf,'w')
    other_ncRNA_gtf_file = open(other_ncRNA_gtf,'w')
    for eachline in GFF_Reader(gtf):
        tr_id = eachline.attr['transcript_id']
        gene_type,exon_length = get_genetype_length(eachline)
        if gene_type in mRNA_tag_list :
            mRNA_gtf_file.write(eachline.get_gff_line())
        elif gene_type in lncRNA_tag_list :
            lncRNA_gtf_file.write(eachline.get_gff_line())
        elif gene_type in lncRNA_candidate_tag_list :
            if tr_id not in candidate_lncRNA_dict :
                candidate_lncRNA_dict.setdefault(tr_id,{})['length'] = exon_length
                candidate_lncRNA_dict.setdefault(tr_id,{})['info'] = [eachline.get_gff_line()]
            else :
                candidate_lncRNA_dict[tr_id]['length'] += exon_length
                candidate_lncRNA_dict[tr_id]['info'].append(eachline.get_gff_line())
        else :
            other_ncRNA_gtf_file.write(eachline.get_gff_line())

    for each_id in candidate_lncRNA_dict :
        if candidate_lncRNA_dict[each_id]['length'] > 200 :
            for eachline in candidate_lncRNA_dict[each_id]['info'] :
                lncRNA_gtf_file.write(eachline)
        else :
            for eachline in candidate_lncRNA_dict[each_id]['info'] :
                other_ncRNA_gtf_file.write(eachline)

    mRNA_gtf_file.close()
    lncRNA_gtf_file.close()
    other_ncRNA_gtf_file.close()
import sys
import os
import argparse
import json
from HTSeq import GFF_Reader

parser = argparse.ArgumentParser()
parser.add_argument('--lnc_classify', help = 'FEELnc_classifier output result', required = True)
parser.add_argument('--lnc_gtf', help = 'lncRNA gtf file', required = True)
parser.add_argument('--output', help = 'output file', required = True)
args = parser.parse_args()

def get_transcript_info(gtf,genename_dict = {}):
    transcript_info_dict = {}
    for eachline in GFF_Reader(gtf):
        gene_id = eachline.attr['gene_id']
        transcript_id  = eachline.attr['transcript_id']
        ## get gene_type info
        if 'gene_type' in eachline.attr :
            gene_type = eachline.attr['gene_type']
        elif 'gene_biotype' in eachline.attr :
            gene_type = eachline.attr['gene_biotype']
        else :
            gene_type = '--'
        ## get gene name info 
        if 'gene_name' in eachline.attr :
            genename = eachline.attr['gene_name']
        elif gene_id in genename_dict :
            genename = genename_dict[gene_id][0]
        else :
            genename = '--'
        ## get location info 
        chrom  = eachline.iv.chrom
        start  = eachline.iv.start + 1
        end    = eachline.iv.end
        strand = eachline.iv.strand
        length = end - start + 1
        if transcript_id in transcript_info_dict :
            if start < transcript_info_dict[transcript_id]['start'] :
                transcript_info_dict[transcript_id]['start'] = start
            if end > transcript_info_dict[transcript_id]['end'] :
                transcript_info_dict[transcript_id]['end'] = end
            transcript_info_dict.setdefault(transcript_id,{})['length'] += length
            transcript_info_dict.setdefault(transcript_id,{})['exon_num'] += 1
        else :            
            transcript_info_dict.setdefault(transcript_id,{})["chrom"] = chrom
            transcript_info_dict.setdefault(transcript_id,{})["start"] = start
            transcript_info_dict.setdefault(transcript_id,{})["end"] = end
            transcript_info_dict.setdefault(transcript_id,{})["strand"] = strand
            transcript_info_dict.setdefault(transcript_id,{})["gene_id"] = gene_id
            transcript_info_dict.setdefault(transcript_id,{})["gene_name"] = genename
            transcript_info_dict.setdefault(transcript_id,{})["gene_description"] = "--"
            transcript_info_dict.setdefault(transcript_id,{})['gene_type'] = gene_type
            transcript_info_dict.setdefault(transcript_id,{})['length'] = length
            transcript_info_dict.setdefault(transcript_id,{})['exon_num'] = 1
            if gene_id in genename_dict :
                transcript_info_dict.setdefault(transcript_id,{})["gene_description"] = genename_dict[gene_id][1]
                gene_description_flag = 1
    return transcript_info_dict    

def get_code(direction,lnc_type,subtype,location):
    if lnc_type == 'intergenic':
        return 'lincRNA'
    else:
        if direction == 'antisense':
            return 'antisense'
        else:
            if location == 'intronic':
                return 'intronic'
            else:
                return 'sense_overlapping'

lnc_info_dict = get_transcript_info(args.lnc_gtf)

lnc_classify_dict = {}
with open(args.lnc_classify) as lnc_classify_info:
    for n,eachline in enumerate(lnc_classify_info):
        eachline_info = eachline.strip().split()
        output_line = '\t'.join(eachline_info[1:])
        if n == 0:
            header = '%s\tlncRNA_class\n' % output_line
        else:
            isBest = eachline_info[0]
            lncRNA_transcript = eachline_info[2]
            direction = eachline_info[5]
            lnc_type = eachline_info[6]
            distance = int(eachline_info[7])
            subtype = eachline_info[8]
            location = eachline_info[9]
            lnc_class = get_code(direction,lnc_type,subtype,location)
            if isBest == '1':
                lnc_classify_dict[lncRNA_transcript] = [lnc_class,output_line]

with open(args.output,'w') as output_info:
    output_info.write(header)
    for each in lnc_info_dict:
        if each in lnc_classify_dict:
            eachline = lnc_classify_dict[each][1]
            lnc_class = lnc_classify_dict[each][0]
            each_out = '%s\t%s\n' % (eachline,lnc_class)
        else:
            gene_id = lnc_info_dict[each]['gene_id']
            other_info = '\t'.join(['--']*7)
            each_out = '%s\t%s\t%s\tlincRNA\n' % (gene_id,each,other_info)
        output_info.write(each_out)





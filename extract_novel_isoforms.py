'''
extract Novel isoforms from assembled gtf
1. class code = 'u,p,i,x'
2. intronic transcript(exon >=2 or antisense)

'''

import sys
import os
from HTSeq import GFF_Reader
import RNAseq_tools

if not len(sys.argv) == 3:
    print '    python ' + sys.argv[0] + ' ref.gtf cuffcompare.merged.gtf > novel.transcript.gtf'
    sys.exit(1)

class_code_dict = {
    'u':'lincRNA',
    'p':'lincRNA',
    'x':'antisense_lncRNA',
    'i':'intronic_lncRNA',
}

ref_gtf = sys.argv[1]
cuffcompare_gtf = sys.argv[2]

ref_tr_dict = RNAseq_tools.get_transcript_info(ref_gtf)
intronic_tr_dict = {}

for eachline in GFF_Reader(cuffcompare_gtf):
    each_tr_id = eachline.attr['transcript_id']
    each_gene_id = eachline.attr['gene_id']
    chrom  = eachline.iv.chrom
    start  = eachline.iv.start + 1
    end    = eachline.iv.end
    strand = eachline.iv.strand
    source = eachline.source
    track_type = eachline.type
    class_code = eachline.attr['class_code']
    if class_code not in class_code_dict:
        continue
    tr_type = class_code_dict[class_code]
    each_track_out = '{chrom}\t{source}\t{track_type}\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{each_gene_id}"; transcript_id "{each_tr_id}"; transcript_type "{tr_type}";'.format(**locals())
    if 'nearest_ref' in eachline.attr:
        nearest_ref_tr = eachline.attr['nearest_ref']
        each_track_out = '{each_track_out} nearest_ref "{nearest_ref_tr}";'.format(**locals())
    if class_code != 'i':
        print each_track_out
    else:
        if each_tr_id not in intronic_tr_dict:
            intronic_tr_dict.setdefault(each_tr_id, {})['exon_number'] = 1
            intronic_tr_dict.setdefault(each_tr_id, {})['strand'] = strand
            intronic_tr_dict.setdefault(each_tr_id, {})['inf'] = [each_track_out]
            intronic_tr_dict.setdefault(each_tr_id, {})['nearest_ref'] = nearest_ref_tr
        else:
            intronic_tr_dict[each_tr_id]['exon_number'] += 1
            intronic_tr_dict[each_tr_id]['inf'].append(each_track_out)

for each_tr_id in intronic_tr_dict:
    each_tr_ref = intronic_tr_dict[each_tr_id]['nearest_ref']
    each_tr_strand = intronic_tr_dict[each_tr_id]['strand']
    each_tr_exon_num = intronic_tr_dict[each_tr_id]['exon_number']
    each_tr_ref_strand = ref_tr_dict[each_tr_ref]['strand']
    if each_tr_exon_num >=2 or each_tr_strand != each_tr_ref_strand:
        for each_inf in intronic_tr_dict[each_tr_id]['inf']:
            print each_inf
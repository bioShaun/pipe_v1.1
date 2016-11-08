import sys
import os
import argparse
from HTSeq import GFF_Reader
import python_tools

mycwd = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument('--gtf',help = 'mRNA lncRNA combined gtf.', required = True)
parser.add_argument('--lncRNA_class',help = 'lncRNA classification file.', required = True)
parser.add_argument('--out_dir', help = 'lncRNA feature file .', required =True)
args = parser.parse_args()

def get_transcript_info(gtf,genename_dict = {}):
    transcript_info_dict = {}
    for eachline in GFF_Reader(gtf):
        if eachline.type != 'exon':
            continue
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
            transcript_info_dict[transcript_id]['exon_start'].append(start)
            transcript_info_dict[transcript_id]['exon_start'].sort()
            transcript_info_dict[transcript_id]['exon_end'].append(end)
            transcript_info_dict[transcript_id]['exon_end'].sort()
            transcript_info_dict[transcript_id]['exon_len'].append(length)          
        else :            
            transcript_info_dict.setdefault(transcript_id,{})["chrom"] = chrom
            transcript_info_dict.setdefault(transcript_id,{})["exon_start"] = [start]
            transcript_info_dict.setdefault(transcript_id,{})["exon_end"] = [end]
            transcript_info_dict.setdefault(transcript_id,{})["strand"] = strand
            transcript_info_dict.setdefault(transcript_id,{})["gene_id"] = gene_id
            transcript_info_dict.setdefault(transcript_id,{})["gene_name"] = genename
            transcript_info_dict.setdefault(transcript_id,{})["gene_description"] = "--"
            transcript_info_dict.setdefault(transcript_id,{})['gene_type'] = gene_type
            transcript_info_dict.setdefault(transcript_id,{})['exon_len'] = [length]
            if gene_id in genename_dict :
                transcript_info_dict.setdefault(transcript_id,{})["gene_description"] = genename_dict[gene_id][1]
                gene_description_flag = 1
    return transcript_info_dict    

def write_list_to_file(obj, fn, header = ''):
    fh = open(fn, 'w')
    if header != '':
        fh.write('%s\n' % header)
    for each_list in obj:
        each_list = [str(each) for each in each_list]
        eachline = '\t'.join(each_list)
        fh.write('%s\n' % eachline)
    fh.close()    

if __name__ == '__main__':
    if not os.path.exists(args.out_dir):
        python_tools.circ_mkdir_unix(args.out_dir)

    tr_dict = get_transcript_info(args.gtf)
    gene_dict = {}
    tr_class_dict = {}
    with open(args.lncRNA_class) as lncRNA_class_info:
        for n, eachline in enumerate(lncRNA_class_info):
            if n != 0:
                eachline_info = eachline.strip().split('\t')
                gene_id = eachline_info[0]
                tr_id = eachline_info[1]
                tr_detail_type = eachline_info[-1]
                tr_type = 'lncRNA'
                tr_class_dict[tr_id] = [tr_type,tr_detail_type]

    lncRNA_exon_len_distribution = os.path.join(args.out_dir,'lncRNA_exon_length.txt')
    mRNA_exon_len_distribution = os.path.join(args.out_dir,'mRNA_exon_length.txt')  
    lncRNA_exon_len_distribution_list = []
    mRNA_exon_len_distribution_list = []

    lncRNA_intron_len_distribution = os.path.join(args.out_dir,'lncRNA_intron_length.txt')
    mRNA_intron_len_distribution = os.path.join(args.out_dir,'mRNA_intron_length.txt')     
    lncRNA_intron_len_distribution_list = []
    mRNA_intron_len_distribution_list = []

    lncRNA_feature = os.path.join(args.out_dir,'lncRNA_feature.txt')
    mRNA_feature = os.path.join(args.out_dir,'mRNA_feature.txt')
    lncRNA_feature_list = []
    mRNA_feature_list = []

    lncRNA_gene_isoform = os.path.join(args.out_dir,'lncRNA_gene_isofrom_number.txt')
    mRNA_gene_isoform = os.path.join(args.out_dir,'mRNA_gene_isofrom_number.txt')
    lncRNA_gene_isoform_list = []
    mRNA_gene_isoform_list = []

    exon_intron_len_distribution_header = 'Chrom\tStart\tEnd\tStrand\tID\tLength\tType'
    tr_feature_info_header = 'Chrom\tStart\tEnd\tStrand\tID\tGene\tLength\tExon_number\tType'
    gene_isoform_info_header = 'ID\tIsoform_number\tIsoform_ID\tType'

    for each_tr in tr_dict:
        each_tr_type = tr_dict[each_tr]['gene_type']
        each_tr_detail_type = each_tr_type
        if each_tr_type == '--':
            if each_tr in tr_class_dict:
                each_tr_type,each_tr_detail_type = tr_class_dict[each_tr]            
            else:
                each_tr_type,each_tr_detail_type = 'TUCP','TUCP'
        each_tr_len = sum(tr_dict[each_tr]['exon_len'])
        each_tr_gene = tr_dict[each_tr]['gene_id']
        each_tr_exon_num = len(tr_dict[each_tr]['exon_len'])
        each_tr_strand = tr_dict[each_tr]['strand']
        each_tr_chrom = tr_dict[each_tr]['chrom']
        each_tr_start = min(tr_dict[each_tr]['exon_start'])
        each_tr_end = max(tr_dict[each_tr]['exon_end'])
        if each_tr_type == 'mRNA' or each_tr_type == 'protein_coding':
            mRNA_feature_list.append([each_tr_chrom, each_tr_start, each_tr_end, each_tr_strand, each_tr, each_tr_gene, each_tr_len, each_tr_exon_num, each_tr_detail_type])
        elif each_tr_type == 'lncRNA':
            lncRNA_feature_list.append([each_tr_chrom, each_tr_start, each_tr_end, each_tr_strand, each_tr, each_tr_gene, each_tr_len, each_tr_exon_num, each_tr_detail_type])            
        else:
            pass
        for n,each_start in enumerate(tr_dict[each_tr]['exon_start']):
            each_end = tr_dict[each_tr]['exon_end'][n]
            if tr_dict[each_tr]['strand'] == '+':
                each_exon_num = n+1
                each_intron_num = n
            else:
                each_exon_num = each_tr_exon_num - n
                each_intron_num = each_tr_exon_num - n
            each_exon_id = '%s.exon.%s' % (each_tr,each_exon_num)
            each_exon_len = each_end - each_start + 1
            if each_tr_type == 'mRNA' or each_tr_type == 'protein_coding':
                mRNA_exon_len_distribution_list.append([each_tr_chrom, each_start, each_end, each_tr_strand, each_exon_id, each_exon_len, each_tr_detail_type])
            elif each_tr_type == 'lncRNA':
                lncRNA_exon_len_distribution_list.append([each_tr_chrom, each_start, each_end, each_tr_strand, each_exon_id, each_exon_len, each_tr_detail_type])                
            else:
                pass
            if n > 0:
                each_intron_start = tr_dict[each_tr]['exon_end'][n-1] + 1
                each_intron_end = each_start - 1 
                each_intron_id = '%s.intron.%s' % (each_tr,each_intron_num)
                each_intron_len = each_intron_end - each_intron_start + 1
                if each_tr_type == 'mRNA' or each_tr_type == 'protein_coding':
                    mRNA_intron_len_distribution_list.append([each_tr_chrom, each_intron_start, each_intron_end, each_tr_strand, each_intron_id, each_intron_len, each_tr_detail_type])
                elif each_tr_type == 'lncRNA':
                    lncRNA_intron_len_distribution_list.append([each_tr_chrom, each_intron_start, each_intron_end, each_tr_strand, each_intron_id, each_intron_len, each_tr_detail_type])
                else:
                    pass
        each_gene_id = tr_dict[each_tr]['gene_id']
        if each_gene_id not in gene_dict:
            gene_dict[each_gene_id] = [[each_tr],[each_tr_type],[each_tr_detail_type]] 
        else:
            gene_dict[each_gene_id][0].append(each_tr)
            gene_dict[each_gene_id][1].append(each_tr_type)
            gene_dict[each_gene_id][2].append(each_tr_detail_type)

for each_gene in gene_dict:
    iso_num = len(gene_dict[each_gene][0])
    iso_info = ','.join(gene_dict[each_gene][0])
    type_list = list(set(gene_dict[each_gene][1]))
    if ('protein_coding' in type_list) or ('mRNA' in type_list):
        each_type = 'protein_coding'
    elif len(type_list) == 1:
        each_type = type_list[0]
    else:
        for each_type in type_list:
            if each_type != 'TUCP' and each_type != 'lncRNA':
                break
        else:
            each_type = 'TUCP'
    if each_type == 'protein_coding' or each_type == 'mRNA':
        mRNA_gene_isoform_list.append([each_gene, iso_num, iso_info, each_type])
    elif each_type == 'lncRNA':
        lncRNA_gene_isoform_list.append([each_gene, iso_num, iso_info, each_type])
    else:
        pass

## sort list
sorted_lncRNA_exon_len_distribution_list = sorted(lncRNA_exon_len_distribution_list, key=lambda x:(x[0],x[1]))
sorted_mRNA_exon_len_distribution_list = sorted(mRNA_exon_len_distribution_list, key=lambda x:(x[0],x[1]))
sorted_lncRNA_intron_len_distribution_list = sorted(lncRNA_intron_len_distribution_list, key=lambda x:(x[0],x[1]))
sorted_mRNA_intron_len_distribution_list = sorted(mRNA_intron_len_distribution_list, key=lambda x:(x[0],x[1]))
sorted_lncRNA_feature_list = sorted(lncRNA_feature_list, key=lambda x:(x[0],x[1]))
sorted_mRNA_feature_list = sorted(mRNA_feature_list, key=lambda x:(x[0],x[1]))
sorted_lncRNA_gene_isoform_list = sorted(lncRNA_gene_isoform_list, key=lambda x:(x[0]))
sorted_mRNA_gene_isoform_list = sorted(mRNA_gene_isoform_list, key=lambda x:(x[0]))

## write list to file
write_list_to_file(sorted_lncRNA_exon_len_distribution_list, lncRNA_exon_len_distribution, exon_intron_len_distribution_header)
write_list_to_file(sorted_mRNA_exon_len_distribution_list, mRNA_exon_len_distribution, exon_intron_len_distribution_header)
write_list_to_file(sorted_lncRNA_intron_len_distribution_list, lncRNA_intron_len_distribution, exon_intron_len_distribution_header)
write_list_to_file(sorted_mRNA_intron_len_distribution_list, mRNA_intron_len_distribution, exon_intron_len_distribution_header)

write_list_to_file(sorted_lncRNA_feature_list, lncRNA_feature, tr_feature_info_header)
write_list_to_file(sorted_mRNA_feature_list, mRNA_feature, tr_feature_info_header)

write_list_to_file(sorted_lncRNA_gene_isoform_list, lncRNA_gene_isoform, gene_isoform_info_header)
write_list_to_file(sorted_mRNA_gene_isoform_list, mRNA_gene_isoform, gene_isoform_info_header)


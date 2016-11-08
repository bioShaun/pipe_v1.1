#! /usr/bin/python

'''

data:2016-3-21

Convert gtf to bed and annotate each feature according to its genomic position (promoter,exon,intron)

'''

import sys
import os
import argparse
from HTSeq import GFF_Reader
import circ

parser = argparse.ArgumentParser()
parser.add_argument('--gtf', help = 'gtf file.', required =True)
parser.add_argument('--bed', help = 'output bed file.', required =True)
args = parser.parse_args()



def sort_list_coordinate(cor_list,rev_flag=False):
    lista = sorted([each[0] for each in cor_list],reverse = rev_flag)
    listb = sorted([each[1] for each in cor_list],reverse = rev_flag)
    return zip(lista,listb)

def exon_to_intron(exon_cor_list):
    exon_num = len(exon_cor_list)
    intron_list = []
    for n in range(1,exon_num) :
        intron_start = exon_cor_list[n-1][1]
        intron_end = exon_cor_list[n][0]
        intron_list.append((intron_start,intron_end))
    return intron_list

def exon_cds_to_utr(exon_cor_list,cds_cor_list,strand):
    if strand == '+' :
        utr5_start = exon_cor_list[0][0]
        utr5_end = cds_cor_list[0][0]
        utr3_start = cds_cor_list[-1][1] + 3
        utr3_end = exon_cor_list[-1][1]
    else :
        utr5_start = cds_cor_list[-1][1]
        utr5_end = exon_cor_list[-1][1]
        utr3_start = exon_cor_list[0][0]
        utr3_end = cds_cor_list[0][0] - 3
    if utr5_start != utr5_end and utr3_start != utr3_end :
        return [(utr5_start,utr5_end)],[(utr3_start,utr3_end)]
    else :
        return None

def get_exon_start_size(exon_cor_list):
    exon_start_list = ['0']
    exon_size_list = []    
    for n,each in enumerate(exon_cor_list) :
        exon_size_list.append(str(each[1]-each[0]))        
        if n == 0 :
            continue          
        else :
            older_start = exon_cor_list[n-1][0]
            exon_start_list.append(str(abs(each[0] - older_start)))
    return ','.join(exon_size_list),','.join(exon_start_list)

def intersect_exon_in_same_tr(exon_cor_list,new_exon):
    for each_exon in exon_cor_list :
        for each_cor in each_exon :
            if each_cor >= new_exon[0] and each_cor <= new_exon[1] :
                return False
    return True

class Transcript_feature:
    def __init__(self):
        self.gtf = None
        self.transcript_dict = {}
        self.tr_pos_map_dict = {}
        self.utr_flag = 0
        self.cds_flag = 0
        self.promoter_dis = 2000
        self.output = None
        self.TYPE = ['exon','intron','cds','five_prime_utr','three_prime_utr']

    def pos_tr_map(self):
        for each_tr in self.transcript_dict :
            start = self.transcript_dict[each_tr]['exon'][0][0]
            chrom = self.transcript_dict[each_tr]['chrom']
            self.tr_pos_map_dict.setdefault(chrom,{})[start] = each_tr

    def get_transcript_dict(self):
        for eachline in GFF_Reader(self.gtf) :
            genomic_featrue = eachline.type.lower()
            if genomic_featrue not in  self.TYPE :
                continue        
            gene_id = eachline.attr['gene_id']
            transcript_id = eachline.attr['transcript_id']
            start = eachline.iv.start
            end = eachline.iv.end
            chrom  = eachline.iv.chrom
            strand = eachline.iv.strand
            self.transcript_dict.setdefault(transcript_id,{})['strand'] = strand
            self.transcript_dict.setdefault(transcript_id,{})['chrom'] = chrom
            self.transcript_dict.setdefault(transcript_id,{})['gene_id'] = gene_id
            if transcript_id in self.transcript_dict and genomic_featrue in self.transcript_dict[transcript_id] :
                if intersect_exon_in_same_tr(self.transcript_dict[transcript_id][genomic_featrue],(start,end)) :
                    self.transcript_dict[transcript_id][genomic_featrue].append((start,end))
            else :
                self.transcript_dict.setdefault(transcript_id,{})[genomic_featrue] = [(start,end)]
            
            genomic_featrue_test = genomic_featrue.upper()            
            if self.cds_flag == 0 and genomic_featrue == 'cds' :
                self.cds_flag = 1
            if self.utr_flag == 0 and 'utr' in genomic_featrue_test :
                self.utr_flag = 1

    def sort_fearture_position(self):
        for each_tr in self.transcript_dict :
            if len(self.transcript_dict[each_tr]['exon']) > 1 :
                self.transcript_dict[each_tr]['exon'] = sort_list_coordinate(self.transcript_dict[each_tr]['exon'])
            if self.cds_flag and 'cds' in self.transcript_dict[each_tr] and len(self.transcript_dict[each_tr]['exon']) > 1:
                self.transcript_dict[each_tr]['cds'] = sort_list_coordinate(self.transcript_dict[each_tr]['cds'])

    def get_intron(self):
        for each_tr in self.transcript_dict :
            if len(self.transcript_dict[each_tr]['exon']) > 1 :
                self.transcript_dict[each_tr]['intron'] = exon_to_intron(self.transcript_dict[each_tr]['exon'])

    def get_promoter(self):
        for each_tr in self.transcript_dict :     
            if self.transcript_dict[each_tr]['strand'] == '+' :
                promoter_start = self.transcript_dict[each_tr]['exon'][0][0] - 2000
                if promoter_start < 0 :
                    promoter_start = 0
                promoter_end = self.transcript_dict[each_tr]['exon'][0][0]
            else :
                promoter_end = self.transcript_dict[each_tr]['exon'][-1][1] + 2000
                promoter_start = self.transcript_dict[each_tr]['exon'][-1][1]
            self.transcript_dict[each_tr]['promoter'] = [(promoter_start,promoter_end)]

    def get_utr(self) :
        if not self.utr_flag and self.cds_flag :
            for each_tr in self.transcript_dict :
                if 'exon' in self.transcript_dict[each_tr] and 'cds' in self.transcript_dict[each_tr] :
                    exon_cor_list = self.transcript_dict[each_tr]['exon']
                    cds_cor_list = self.transcript_dict[each_tr]['cds']
                    strand = self.transcript_dict[each_tr]['strand']
                    if exon_cds_to_utr(exon_cor_list,cds_cor_list,strand) :
                        self.transcript_dict[each_tr]['five_prime_utr'],self.transcript_dict[each_tr]['three_prime_utr'] = exon_cds_to_utr(exon_cor_list,cds_cor_list,strand)

    def output_bed_info(self,tr_id,cor_list,genomic_featrue):
        feature_num = 0
        chrom = self.transcript_dict[tr_id]['chrom']
        strand = self.transcript_dict[tr_id]['strand']
        gene_id = self.transcript_dict[tr_id]['gene_id']
        feature_info_list = []    
        for each_feature in cor_list :
            feature_num += 1
            feature_start = each_feature[0]
            feature_end = each_feature[1]
            feature_length = feature_end - feature_start
            if genomic_featrue in ['exon','intron','cds'] :
                feature_info_list.append('{chrom}\t{feature_start}\t{feature_end}\t{tr_id}|{genomic_featrue}.{feature_num}\t0\t{strand}\t{feature_start}\t{feature_end}\t0\t1\t{feature_length},\t0,\t{gene_id}'.format(**locals()))                
            else :
                feature_info_list.append('{chrom}\t{feature_start}\t{feature_end}\t{tr_id}|{genomic_featrue}\t0\t{strand}\t{feature_start}\t{feature_end}\t0\t1\t{feature_length},\t0,\t{gene_id}'.format(**locals()))                

        return feature_info_list

    def output_bed(self):
        self.get_transcript_dict()
        self.sort_fearture_position()
        self.get_intron()
        self.get_promoter()
        self.get_utr()
        self.pos_tr_map()
        feature_info_list = []
        chrom_list = sorted(self.tr_pos_map_dict.keys())
        for each_chr in chrom_list :
            sorted_pos_list = sorted(self.tr_pos_map_dict[each_chr].keys())
            for each_pos in sorted_pos_list :
                tr_id = self.tr_pos_map_dict[each_chr][each_pos]
                exon_cor_list = self.transcript_dict[tr_id]['exon'][:]
                gene_id = self.transcript_dict[tr_id]['gene_id']
                strand = self.transcript_dict[tr_id]['strand']
                exon_num = len(exon_cor_list)
                tr_start = exon_cor_list[0][0]
                tr_end = exon_cor_list[-1][1]
                if strand == '-' :
                    exon_cor_list = sort_list_coordinate(exon_cor_list,True)
                tr_exon_size,tr_exon_length = get_exon_start_size(exon_cor_list)
                ## transcript info
                feature_info_list.append('{each_chr}\t{tr_start}\t{tr_end}\t{tr_id}|transcript\t0\t{strand}\t{tr_start}\t{tr_end}\t0\t{exon_num}\t{tr_exon_size},\t{tr_exon_length}\t{gene_id}'.format(**locals()))                
                ## output_promoter
                promoter_cor_list = self.transcript_dict[tr_id]['promoter'][:]
                promoter_info = self.output_bed_info(tr_id,promoter_cor_list,'promoter')        
                feature_info_list.extend(promoter_info)                
                ## output_exon                                
                each_exon_info_list = self.output_bed_info(tr_id,exon_cor_list,'exon')
                feature_info_list.extend(each_exon_info_list)
                ## output_cds 
                if 'cds' in self.transcript_dict[tr_id] :
                    cds_cor_list = self.transcript_dict[tr_id]['cds'][:]
                    if strand == '-' :
                        cds_cor_list = sort_list_coordinate(cds_cor_list,True)
                    each_cds_info_list = self.output_bed_info(tr_id,cds_cor_list,'cds')
                    feature_info_list.extend(each_cds_info_list)
                ## output intron
                if 'intron' in self.transcript_dict[tr_id] :
                    intron_cor_list = self.transcript_dict[tr_id]['intron']
                    if strand == '-' :
                        intron_cor_list = sort_list_coordinate(intron_cor_list,True)
                    each_intron_info_list = self.output_bed_info(tr_id,intron_cor_list,'intron')
                    feature_info_list.extend(each_intron_info_list)                    
                ## output utr
                if 'five_prime_utr' in self.transcript_dict[tr_id] :
                    utr5_cor_list = self.transcript_dict[tr_id]['five_prime_utr'][:]                    
                    utr5_info = self.output_bed_info(tr_id,utr5_cor_list,'five_prime_utr')                            
                    feature_info_list.extend(utr5_info)               
                if 'three_prime_utr' in self.transcript_dict[tr_id] :
                    utr3_cor_list = self.transcript_dict[tr_id]['three_prime_utr'][:]
                    utr3_info = self.output_bed_info(tr_id,utr3_cor_list,'three_prime_utr')    
                    feature_info_list.extend(utr3_info)
        circ.write_obj_to_file(feature_info_list,self.output)

def main():
    my_Transcript_feature = Transcript_feature()
    my_Transcript_feature.gtf = args.gtf
    my_Transcript_feature.output = args.bed
    my_Transcript_feature.output_bed()
#     my_Transcript_feature.get_transcript_dict()

if __name__ == '__main__':
    main()

# if __name__ == '__main__':
#     my_Transcript_feature = Transcript_feature()
#     my_Transcript_feature.gtf = args.gtf
#     my_Transcript_feature.output = args.bed
#     #my_Transcript_feature.output_bed()
#     my_Transcript_feature.get_transcript_dict()
#     my_Transcript_feature.sort_fearture_position()
#     my_Transcript_feature.pos_tr_map()
#     my_Transcript_feature.get_intron()
#     my_Transcript_feature.get_promoter()
#     my_Transcript_feature.get_utr()

#     feature_info_list = []
#     chrom_list = sorted(my_Transcript_feature.tr_pos_map_dict.keys())
#     for each_chr in chrom_list :
#         sorted_pos_list = sorted(my_Transcript_feature.tr_pos_map_dict[each_chr].keys())
#         for each_pos in sorted_pos_list :
#             tr_id = my_Transcript_feature.tr_pos_map_dict[each_chr][each_pos]
#             exon_cor_list = my_Transcript_feature.transcript_dict[tr_id]['exon'][:]
#             gene_id = my_Transcript_feature.transcript_dict[tr_id]['gene_id']
#             strand = my_Transcript_feature.transcript_dict[tr_id]['strand']
#             exon_num = len(exon_cor_list)
#             tr_start = exon_cor_list[0][0]
#             tr_end = exon_cor_list[-1][1]
#             if strand == '-' :
#                 exon_cor_list = sort_list_coordinate(exon_cor_list,True)
#             tr_exon_size,tr_exon_length = get_exon_start_size(exon_cor_list)
#             ## transcript info
#             feature_info_list.append('{each_chr}\t{tr_start}\t{tr_end}\t{tr_id}|transcript\t0\t{strand}\t{tr_start}\t{tr_end}\t0\t{exon_num}\t{tr_exon_size},\t{tr_exon_length}\t{gene_id}'.format(**locals()))                
#             ## output_promoter
#             promoter_cor_list = my_Transcript_feature.transcript_dict[tr_id]['promoter'][:]
#             promoter_info = my_Transcript_feature.output_bed_info(tr_id,promoter_cor_list,'promoter')        
#             feature_info_list.extend(promoter_info)                
#             ## output_exon                                
#             exon_num = 0
#             each_exon_info_list = my_Transcript_feature.output_bed_info(tr_id,exon_cor_list,'exon')
#             feature_info_list.extend(each_exon_info_list)
#             ## output_cds 
#             if 'cds' in my_Transcript_feature.transcript_dict[tr_id] :
#                 cds_cor_list = my_Transcript_feature.transcript_dict[tr_id]['cds'][:]
#                 if strand == '-' :
#                     cds_cor_list = sort_list_coordinate(cds_cor_list,True)
#                 each_cds_info_list = my_Transcript_feature.output_bed_info(tr_id,cds_cor_list,'cds')
#                 feature_info_list.extend(each_cds_info_list)
#             ## output intron
#             if 'intron' in my_Transcript_feature.transcript_dict[tr_id] :
#                 intron_cor_list = my_Transcript_feature.transcript_dict[tr_id]['intron']
#                 if strand == '-' :
#                     intron_cor_list = sort_list_coordinate(intron_cor_list,True)
#                 each_intron_info_list = my_Transcript_feature.output_bed_info(tr_id,intron_cor_list,'intron')
#                 feature_info_list.extend(each_intron_info_list)                    
#             ## output utr
#             if 'five_prime_utr' in my_Transcript_feature.transcript_dict[tr_id] :
#                 utr5_cor_list = my_Transcript_feature.transcript_dict[tr_id]['five_prime_utr'][:]
#                 utr5_info = my_Transcript_feature.output_bed_info(tr_id,utr5_cor_list,'five_prime_utr')                  
#                 feature_info_list.extend(utr5_info)                
#             if 'three_prime_utr' in my_Transcript_feature.transcript_dict[tr_id] :
#                 utr3_cor_list = my_Transcript_feature.transcript_dict[tr_id]['three_prime_utr'][:]
#                 utr3_info = my_Transcript_feature.output_bed_info(tr_id,utr3_cor_list,'three_prime_utr')
#                 feature_info_list.extend(utr3_info)
#     circ.write_obj_to_file(feature_info_list,my_Transcript_feature.output)


#! /usr/bin/python

#coding=utf-8

import sys
import HTSeq
import argparse

if len(sys.argv) != 3 :
    print "python "+sys.argv[0]+" gff gtf"
    sys.exit(0)

genes_dic = {}
transcirpts_dic = {}
cds_gtf_dic = {}
exon_gtf_dic = {}

for eachline in HTSeq.GFF_Reader(sys.argv[1]):
    if eachline.type == "gene" :
        if eachline.attr['ID'] not in genes_dic :
            gene_id = eachline.attr['ID']
            gene_name = "--"
            genes_dic[eachline.attr['ID']] = []
            genes_dic[eachline.attr['ID']].append(gene_id)
            genes_dic[eachline.attr['ID']].append(gene_name)
            if eachline.attr.has_key('gene') :
                genes_dic[eachline.attr['ID']][1] = eachline.attr['gene']
            elif eachline.attr.has_key('Name') :
                genes_dic[eachline.attr['ID']][1] = eachline.attr['Name']
            if eachline.attr.has_key('locus_tag') :
                genes_dic[eachline.attr['ID']][0] = eachline.attr['locus_tag']

    elif eachline.type == "mRNA" :
        if eachline.attr['ID'] not in transcirpts_dic and (eachline.attr.has_key('Parent')):
            transcirpts_dic[eachline.attr['ID']] = eachline.attr['Parent']
        elif eachline.attr['ID'] not in transcirpts_dic and (not eachline.attr.has_key('Parent')):
            transcirpts_dic[eachline.attr['ID']] = 1

    elif eachline.type == "CDS" :
        cds_parent = eachline.attr['Parent']
        gene_id,transcript_id = cds_parent,cds_parent
        gene_type = "mRNA"
        gene_name = "--"
        if transcript_id in transcirpts_dic :
            gene_id = transcirpts_dic[transcript_id]
            if gene_id != 1 and gene_id in genes_dic:
                gene_name = genes_dic[gene_id][1]
                gene_id = genes_dic[gene_id][0]
            else :
                gene_id = transcript_id
            output_gene_id,output_transcript_id = gene_id,transcript_id

        elif gene_id in genes_dic :
            gene_name = genes_dic[gene_id][1]
            output_gene_id,output_transcript_id = genes_dic[gene_id][0],genes_dic[gene_id][0]

        all_cds_info = eachline.get_gff_line().strip().split("\t")
        cds_info_output = all_cds_info[0:8]
        cds_info_output[2] = "exon"
        cds_attribute = "gene_id \"{output_gene_id}\"; transcript_id \"{output_transcript_id}\"; gene_name \"{gene_name}\"; gene_type \"{gene_type}\";\n".format(**locals())
        cds_info_output.append(cds_attribute)
        cds_info_output_line = "\t".join(cds_info_output)
        if transcript_id not in cds_gtf_dic :
            cds_gtf_dic[transcript_id] = []
            cds_gtf_dic[transcript_id].append(cds_info_output_line)
        else :
            cds_gtf_dic[transcript_id].append(cds_info_output_line)

    elif eachline.type == "exon" :
        exon_parent = eachline.attr['Parent']
        gene_id,transcript_id = exon_parent,exon_parent
        if eachline.attr.has_key('gbkey') :
            gene_type = eachline.attr['gbkey']
        else :
            gene_type = "mRNA"
        gene_name = "--"
        if transcript_id in transcirpts_dic :
            gene_id = transcirpts_dic[transcript_id]
            if gene_id != 1 and gene_id in genes_dic:
                gene_name = genes_dic[gene_id][1]
                gene_id = genes_dic[gene_id][0]
            else :
                gene_id = transcript_id
            output_gene_id,output_transcript_id = gene_id,transcript_id

        elif gene_id in genes_dic :
            gene_name = genes_dic[gene_id][1]
            output_gene_id,output_transcript_id = genes_dic[gene_id][0],genes_dic[gene_id][0]

        all_exon_info = eachline.get_gff_line().strip().split("\t")
        exon_info_output = all_exon_info[0:8]
        exon_attribute = "gene_id \"{output_gene_id}\"; transcript_id \"{output_transcript_id}\"; gene_name \"{gene_name}\"; gene_type \"{gene_type}\";\n".format(**locals())
        exon_info_output.append(exon_attribute)
        exon_info_output_line = "\t".join(exon_info_output)
        if transcript_id not in exon_gtf_dic :
            exon_gtf_dic[transcript_id] = []
            exon_gtf_dic[transcript_id].append(exon_info_output_line)
        else :
            exon_gtf_dic[transcript_id].append(exon_info_output_line)

all_id = genes_dic.keys()
all_id.extend(transcirpts_dic.keys())
all_id = list(set(all_id))

output = open(sys.argv[2],"w")
for each_transcirpt in all_id :
    if each_transcirpt in exon_gtf_dic:
        for each_exon in exon_gtf_dic[each_transcirpt] :
            output.write(each_exon)
    elif each_transcirpt in cds_gtf_dic :
        for each_cds in cds_gtf_dic[each_transcirpt] :
            output.write(each_cds)

output.close()

'''
data:2016-3-16

lncRNA filtering

step1 --- length filter (> 200bp)
step2 --- protein coding transcripts filter

'''
import sys
import os
import argparse
from HTSeq import GFF_Reader
import json
import glob
import python_tools

parser = argparse.ArgumentParser()
parser.add_argument('--assembly', help = 'assembled GTF file', required = True)
parser.add_argument('--ref', help = 'ref GTF file', required = True)
parser.add_argument('--out_dir', help = 'Output directory.', required = True)
parser.add_argument('--quant_file', help = 'transcripts expression tables, seprated with "," ',required = True)
parser.add_argument('--replicate_info', help = 'file with 1st column sample id, 2nd column replicate id', default = '')
args = parser.parse_args()

CUFF_BIN = "/usr/bin/"
GTF2BED = '/home/lxgui/scripts/gtf2Bed.pl'
GTF_CLASSIFY = '/home/lxgui/scripts/gtf_classification.py'
LNC_LENGTH = 200
OVERLAPPING_FLAGS = '=,c,j,e,o'.split(',')
LNCRNA_TYPE = ('lincRNA','antisense','antisense_intronic','sense_intronic','ncRNA_host')

def mkdir_unix(path):
    cmd = 'mkdir -p %s' % path
    if not os.path.isdir(path):
        os.system(cmd)

def get_class_code(gffline) :
    old_tr_id = gffline.attr['oId']
    class_code = gffline.attr['class_code']
    nearest_ref = ''
    if 'nearest_ref' in gffline.attr :
        nearest_ref = gffline.attr['nearest_ref']
    return old_tr_id,class_code,nearest_ref

def gtf2bed(gtf,name,out_dir):  
    gtf2bed_perl = GTF2BED
    sorted_gtf = os.path.join(out_dir,'tmp.gtf')
    bedfile = os.path.join(out_dir,'tmp.bed')
    sorted_bed = os.path.join(out_dir,'%s.bed' % name)
    # os.system('sort -k1,1 -k4,4n {gtf} > {sorted_gtf}'.format(**locals()))
    # os.system('perl {gtf2bed_perl} {sorted_gtf} > {bedfile}'.format(**locals()))
    # os.system('sort -k1,1 -k2,2n {bedfile} > {sorted_bed}'.format(**locals()))
    # os.system('rm {bedfile} {sorted_gtf}'.format(**locals()))
    return sorted_bed

def get_basic_info(bedfile):
    tr_dict = {}
    with open(bedfile,'r') as bedfile_info :
        for eachline in bedfile_info :
            if eachline.strip() != "" :
                each_tr_info = eachline.strip().split('\t')
                tr_id = each_tr_info[3]
                tr_dict.setdefault(tr_id,{})['chrom'] = each_tr_info[0]
                tr_dict.setdefault(tr_id,{})['start'] = int(each_tr_info[1])
                tr_dict.setdefault(tr_id,{})['end'] = int(each_tr_info[2])
                tr_dict.setdefault(tr_id,{})['strand'] = each_tr_info[5]
                tr_dict.setdefault(tr_id,{})['exon_num'] = int(each_tr_info[9])
                tr_dict.setdefault(tr_id,{})['tr_length'] = sum([ int(each) for each in each_tr_info[10].split(',') if each != ''])
    return tr_dict

def ref_classification(gtf,out_dir):
    ref_dir = os.path.join(out_dir,'ref')
    mkdir_unix(ref_dir)
    gtf_classify = GTF_CLASSIFY
    # os.system('python {gtf_classify} --gtf {gtf} --out_dir {ref_dir}'.format(**locals()))
    return ref_dir

def compare_with_ref(ref,qurey,out_dir):
    combined_gtf = ''
    if os.path.getsize(ref) :
        cuffcompare = os.path.join(CUFF_BIN,'cuffcompare')
        name = os.path.basename(ref)
        # os.system('{cuffcompare} -r {ref} -R {qurey} -T -o {out_dir}/{name}'.format(**locals()))
        combined_gtf = '{out_dir}/{name}.combined.gtf'.format(**locals())
    return combined_gtf

def store_into_json(obj,fn):
    with open(fn,'w') as file_info :
        json.dump(obj,file_info)

def load_from_json(fn) :
    with open(fn,'r') as file_info :
        obj = json.load(file_info) 
    return obj

def sample_to_replicate(replicate_info):
    sample_to_replicate_dict = {}
    with open(replicate_info,'r') as replicate_info_file :
        for eachline in replicate_info_file :
            eachline_info = eachline.strip().split('\t')
            rep_id = eachline_info[1]
            sample_id = eachline_info[0]
            if sample_id not in sample_to_replicate_dict :
                sample_to_replicate_dict[sample_id] = rep_id
            else :
                sys.exit('duplicate sample id')  
    return sample_to_replicate_dict             

class lncRNA_filter :
    def __init__(self):
        self.assemlby_dict = {}
        self.assembly_dict_json = None
        self.assembly_gtf = None
        self.ref_dict = None
        self.out_dir = None
        self.combined_gtf_list = []
        self.length = LNC_LENGTH
        self.overlap_flags = OVERLAPPING_FLAGS
        self.lncRNA_type = LNCRNA_TYPE
        self.novel_transcript_gtf = None
        self.quant_file_list = []
        self.transcript_exp_dict = {}
        self.exp_novel_transcript_gtf = None

    def add_compare_info(self,combined_gtf):
        if 'mRNA' in combined_gtf :
            for eachline in GFF_Reader(combined_gtf):
                tr_id,class_code,nearest_ref = get_class_code(eachline)
                self.assemlby_dict[tr_id]['nearest_ref'] = nearest_ref
                if class_code in ['u','p'] :
                    self.assemlby_dict[tr_id]['status'] = 'lincRNA'
                elif class_code == 'x' :
                    self.assemlby_dict[tr_id]['status'] = 'antisense'                    
                elif class_code == 'i' :
                    if self.assemlby_dict[tr_id]['strand'] != self.ref_dict[nearest_ref]['strand'] :
                        self.assemlby_dict[tr_id]['status'] = 'antisense_intronic'
                    else :
                        if self.assemlby_dict[tr_id]['exon_num'] > 1 :
                            self.assemlby_dict[tr_id]['status'] = 'sense_intronic'
                        else :
                            self.assemlby_dict[tr_id]['status'] = 'backgroud'
                else :
                    self.assemlby_dict[tr_id]['status'] = 'protein_coding'

        elif 'other_ncRNA' in combined_gtf :
            for eachline in GFF_Reader(combined_gtf):
                tr_id,class_code,nearest_ref = get_class_code(eachline)
                if class_code in self.overlap_flags and 'status' in self.assemlby_dict[tr_id]  :
                    self.assemlby_dict[tr_id]['status'] = 'ncRNA_host'

        elif 'lncRNA' in combined_gtf :
            for eachline in GFF_Reader(combined_gtf):
                tr_id,class_code,nearest_ref = get_class_code(eachline)
                if class_code in self.overlap_flags and 'status' in self.assemlby_dict[tr_id] :
                    self.assemlby_dict[tr_id]['status'] = 'Annotated_lncRNA'                  

    def get_novel_transcript(self):
        ## annotate mRNA related transcript
        for each_gtf in self.combined_gtf_list :
            self.add_compare_info(each_gtf)
        ## filter transcripts
        self.assembly_dict_json = os.path.join(self.out_dir,'assembly_tr_info.json')
        store_into_json(self.assemlby_dict,self.assembly_dict_json)
        ## output filtered transcript gtf
        self.novel_transcript_gtf = os.path.join(self.out_dir,'novel_transcript.gtf')
        output = open(self.novel_transcript_gtf,'w')
        for eachline in GFF_Reader(self.assembly_gtf):
            tr_id = eachline.attr['transcript_id']
            ## filter length
            if self.assemlby_dict[tr_id]['tr_length'] <= self.length :
                continue
            ## annotation filter
            elif 'status' not in self.assemlby_dict[tr_id] : 
                continue
            elif self.assemlby_dict[tr_id]['status'] not in self.lncRNA_type :
                continue
            else :
                output.write(eachline.get_gff_line())
        output.close()

    def store_exp(self,quant_file,sample_to_replicate_dict) :
        with open(quant_file,'r') as quant_file_info :
            for n,eachline in enumerate(quant_file_info) :
                eachline_info = eachline.strip().split('\t')                
                if n == 0 :
                    sample_id_list = eachline_info[1:]
                else :
                    tr_id = eachline_info[0]
                    exp_list = eachline_info[1:]
                    for n in range(len(exp_list)) :
                        sample_id = '_'.join(sample_id_list[n].split('_')[:-1])
                        sample_exp = float(exp_list[n])
                        replicate = sample_to_replicate_dict[sample_id]
                        if tr_id not in self.transcript_exp_dict :
                            self.transcript_exp_dict.setdefault(tr_id,{})[replicate] = [sample_exp]
                        else :
                            python_tools.add_dict_value(self.transcript_exp_dict[tr_id],replicate,sample_exp)
       
    def get_expressed_transcript(self):
        self.novel_transcript_gtf = os.path.join(self.out_dir,'novel_transcript.gtf')
        sample_to_replicate_dict = sample_to_replicate(args.replicate_info)
        for each_quant_file in self.quant_file_list :
            self.store_exp(each_quant_file,sample_to_replicate_dict)
        self.exp_novel_transcript_gtf = os.path.join(self.out_dir,'exp_novel_transcript.gtf')
        output = open(self.exp_novel_transcript_gtf,'w')
        exp_flag_dict = {}
        for eachline in GFF_Reader(self.novel_transcript_gtf) :
            tr_id = eachline.attr['transcript_id']
            if tr_id in exp_flag_dict and exp_flag_dict[tr_id] :
                output.write(eachline.get_gff_line())
                continue    
            elif tr_id not in self.transcript_exp_dict :
                continue
            for each_rep in self.transcript_exp_dict[tr_id] :
                exp_flag_dict[tr_id] = False   
                tr_exp = min(self.transcript_exp_dict[tr_id][each_rep])
                # tr_exp = python_tools.Median(self.transcript_exp_dict[tr_id][each_rep])
                tr_exon_num = self.assemlby_dict[tr_id]['exon_num']
                if tr_exon_num == 1 and tr_exp >= 2 :
                    output.write(eachline.get_gff_line())                 
                    break
                elif tr_exp >= 0.5 :
                    output.write(eachline.get_gff_line())
                    break
            else :
                exp_flag_dict[tr_id] = True

        output.close()

if __name__ == '__main__':
    mkdir_unix(args.out_dir)
    ##-- get bed file --##
    assembly_bed = gtf2bed(args.assembly,'assembly',args.out_dir)
    ref_bed = gtf2bed(args.ref,'ref',args.out_dir)
    ##-- get assembly and ref transcript infomation --##
    assembly_tr_dict = get_basic_info(assembly_bed)
    ref_tr_dict = get_basic_info(ref_bed)
    ##-- classify ref gtf to mRNA,lncRNA and other ncRNA --##
    ref_gtf_dir = ref_classification(args.ref,args.out_dir)
    ##-- compare assembly gtf with ref --##
    combined_gtf_dir = os.path.join(args.out_dir,'compare_dir')    
    mkdir_unix(combined_gtf_dir)
    ref_gtf_list = [os.path.join(ref_gtf_dir,each) for each in os.listdir(ref_gtf_dir)]
    combined_gtf_list = []
    for each_gtf in ref_gtf_list :
        combined_gtf_file = compare_with_ref(each_gtf,args.assembly,combined_gtf_dir)
        if combined_gtf_file :
            combined_gtf_list.append(combined_gtf_file)
    ##-- filtering assembly gtf --##
    my_lncRNA_filter = lncRNA_filter()
    # my_lncRNA_filter.assemlby_dict = assembly_tr_dict
    my_lncRNA_filter.assembly_gtf = args.assembly
    my_lncRNA_filter.ref_dict = ref_tr_dict
    my_lncRNA_filter.out_dir = args.out_dir
    my_lncRNA_filter.combined_gtf_list = combined_gtf_list
    ##-- get novel transcripts --##
    # my_lncRNA_filter.get_novel_transcript()
    ##-- get expressed novel transcripts --##  
    my_lncRNA_filter.quant_file_list = args.quant_file.split(',')
    assembled_json = os.path.join(args.out_dir,"assembly_tr_info.json")
    my_lncRNA_filter.assemlby_dict = load_from_json(assembled_json)
    sample_to_replicate_dict = sample_to_replicate(args.replicate_info)
    my_lncRNA_filter.get_expressed_transcript()


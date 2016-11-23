#coding=utf-8

'''
data:2016-3-26

RNAseq tools and functions

function list:

last modification:2016-3-26
'''

import sys
import os
import python_tools
import gzip
from HTSeq import GFF_Reader
import hashlib   
import logging
import gzip
import glob
import re
import json
import ConfigParser

## database
KO_PEP_DIR = '/home/lxgui/test_kobas/seq_pep/'

## softwares and scripts absolute path
R_BIN = '/usr/local/bin/'
R_BIN_3_1 = '/home/public/R_packages/R-3.1.0/bin/'
BLAST_BIN = '/usr/bin/'
KALLISTO = '/home/public/software/kallisto-0.42.4/kallisto'
TH_CORE_NUM = 24
LC_228_CORE_NUM = 64
LC_34_CORE_NUM = 4
SERVER_CORE_DICT = {
    'server228':24,
    'server34':4,
    'th':24,
}


GOSEQ = '/home/lxgui/scripts/goseq.R'
TOPGO = '/home/lxgui/scripts/topGO_plot.R'
# ENRICH_BAR = '/home/lxgui/scripts/enrich_bar.2016-3-31.R'

TREAT_KEGG_TABLE = '/home/lxgui/scripts/treat_KEGG_table.py'
EXTRACT_INF_BY_ID = '/home/lxgui/scripts/extract_info_by_id.py'

ENRICH_BAR = '/home/lxgui/scripts/enrich_bar.2016-4-13_n.R'
PATHVIEW = '/home/lxgui/scripts/kegg_pathview.py'
PATHVIEW_SRNA = '/home/lxgui/scripts/kegg_pathview_sRNA.py'
PATHVIEW_CK = '/home/lxgui/scripts/check_kegg_pathway.py'

READS_QUALITY_PLOT = '/home/lxgui/scripts/qc/plot/reads_quality_barplot.R'
EXP_PLOT = '/home/lxgui/scripts/quant/plot/expression_boxplot_v1.2.R'
DIFF_TABLE_TREAT = '/home/lxgui/scripts/diff_table_add_header.py'
VOLCANO_PLOT = '/home/lxgui/scripts/Volcano_Plot_20160406.R'
PCA_PLOT = '/home/lxgui/scripts/quant/plot/PCA.R'
CLUSTER_LINE = '/home/lxgui/scripts/quant/plot/cluster_plot.R'

GO_RESULT_ANNO = '/home/lxgui/scripts/GO_add_diff_gene.py'
QUANT_ANNO = '/home/lxgui/scripts/add_gene_anno_v2.py'
FASTQC_SUMMERY = '/home/lxgui/scripts/fastqc_get_data_info_pipe.py'

def get_transcript_info(gtf,genename_dict = {}):
    transcript_info_dict = {}
    for eachline in GFF_Reader(gtf):
        gene_id = eachline.attr['gene_id']
        if 'transcript_id' not in eachline.attr:
            continue
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

def get_gene_info(transcript_info_dict) :
    gene_info_dict = {}
    for each_tr in transcript_info_dict :
        gene_id = transcript_info_dict[each_tr]["gene_id"]
        each_tr_len = transcript_info_dict[each_tr]['length']
        each_tr_exon = transcript_info_dict[each_tr]['exon_num']
        each_tr_tss = transcript_info_dict[each_tr]['start']
        each_tr_tts = transcript_info_dict[each_tr]['end']
        each_gene_strand = transcript_info_dict[each_tr]['strand']
        each_gene_chrom = transcript_info_dict[each_tr]['chrom']
        if gene_id in gene_info_dict :
            gene_info_dict[gene_id]['transcript_id'].append(each_tr)
            gene_info_dict[gene_id]['transcript_len'].append(each_tr_len)
            gene_info_dict[gene_id]['tss'].append(each_tr_tss)
            gene_info_dict[gene_id]['tts'].append(each_tr_tts)
            gene_info_dict[gene_id]['exon_num'].append(each_tr_exon)
        else :
            gene_info_dict.setdefault(gene_id,{})['transcript_id'] = [each_tr]
            gene_info_dict.setdefault(gene_id,{})['transcript_len'] = [each_tr_len]
            gene_info_dict.setdefault(gene_id,{})['tss'] = [each_tr_tss]
            gene_info_dict.setdefault(gene_id,{})['tts'] = [each_tr_tts]
            gene_info_dict.setdefault(gene_id,{})['exon_num'] = [each_tr_exon]
            gene_info_dict.setdefault(gene_id,{})['strand'] = each_gene_strand
            gene_info_dict.setdefault(gene_id,{})['chrom'] = each_gene_chrom
    return gene_info_dict

def enrich_bar_plot(compare_name,output,outname,out_dir) :
    cmd = '%s/Rscript %s %s %s %s %s' % (R_BIN_3_1,ENRICH_BAR,compare_name,output,outname,out_dir)
    return cmd

def get_fq_sample_name(fqname) :
    if '_combined_R1.fastq.gz' in fqname :
        return 'R1',fqname.split('_combined_R1.fastq.gz')[0]
    else :
        return 'R2',fqname.split('_combined_R2.fastq.gz')[0]

def md5sum(fname):
    """ 
    计算文件的MD5值
    
    """
    def read_chunks(fh):
        fh.seek(0)
        chunk = fh.read(8096)
        while chunk:
            yield chunk
            chunk = fh.read(8096)
        else: #最后要将游标放回文件开头
            fh.seek(0)
    m = hashlib.md5()
    if isinstance(fname, basestring) \
        and os.path.exists(fname):
        with open(fname, "rb") as fh:
            for chunk in read_chunks(fh):
                m.update(chunk)
    #上传的文件缓存 或 已打开的文件流
    elif fname.__class__.__name__ in ["StringIO", "StringO"] \
        or isinstance(fname, file):
        for chunk in read_chunks(fname):
            m.update(chunk)
    else:
        return ""
    return m.hexdigest()

class treat_fq :
    def __init__(self):
        self.raw_dir = None
        self.clean_dir = None
        self.mapfile = None
        self.fq_flag = 'PE'
        self.fq_dict = {}
        self.lib_sampl_dict = {}
        self.sample_name_order = []

    def lib_name_to_sample_name(self):
        if self.mapfile :
            self.lib_sampl_dict = python_tools.table_to_dict(self.mapfile,2,1,'\t',False)

    def get_sample_name_order(self) :
        if self.mapfile :
            self.sample_name_order = [each.strip().split('\t')[0] for each in open(self.mapfile,'r')]

    def get_sample_name(self,fq_file):
        fq_name = os.path.basename(fq_file)
        if 'fastq' in fq_name :
            sample_name_tmp = fq_name.split('.fastq')[0]
        elif 'fq' in fq_name :
            sample_name_tmp = fq_name.split('.fq')[0]
        else :
            sys.exit('"fq" and "fastq" not in fq file name ! : ' % fq_file)
        if self.fq_flag == 'SE' :
            sample_name = sample_name_tmp
            if self.lib_sampl_dict :
                sample_name = self.lib_sampl_dict[sample_name]
            self.fq_dict[sample_name] = fq_file   
        else :
            sample_name = '_'.join(sample_name_tmp.split('_')[0:-1])
            fq_order = sample_name_tmp.split('_')[-1]
            if self.lib_sampl_dict :
                sample_name = self.lib_sampl_dict[sample_name]
            self.fq_dict.setdefault(sample_name,{})[fq_order] = fq_file
        return sample_name    

    def get_fq_len(self,fq_file) :
        if fq_file.endswith('.gz') :
            fq_file_info = gzip.open(fq_file, 'rb')
        else :
            fq_file_info = open(fq_file, 'r')
        fq_file_info.readline()
        test_fq_len = len(fq_file_info.readline().strip())
        fq_file_info.close()
        return test_fq_len

    def raw_to_clean(self) :
        pass

    def store_clean_fq(self):
        fq_file_list = os.listdir(self.clean_dir)
        for each_fq in fq_file_list :
            each_fq_file = os.path.join(self.clean_dir,each_fq)
            self.get_sample_name(each_fq_file)
        
class run_kallisto :
    def __init__(self):
        self.kallisto = KALLISTO
        self.transcript_fa = None
        self.index = None
        self.fq_list = []
        self.out_dir = None
        self.fq_length = None
        self.sd = 20
        self.thread = 4
        self.sample_list = []

    def kallisto_index(self):
        cmd = '%s index -i %s %s' % (self.kallisto,self.index,self.transcript_fa)
        python_tools.circ_call_process(cmd)
        run_log = os.path.join(self.out_dir,'index.cmd.log')
        python_tools.write_obj_to_file(cmd,run_log)
        return cmd

    def kallisto_quant(self):
        fq_cmd = ' '.join(self.fq_list)
        if len(self.fq_list) > 1 :  
            cmd = '%s quant -i %s -o %s %s' % (self.kallisto,self.index,self.out_dir,fq_cmd)          
            python_tools.circ_call_process(cmd)
        else :
            single_length = 2*self.fq_length
            cmd = '%s quant -i %s -o %s --single -l %s -s %s --plaintext -t %s %s ' % (self.kallisto,self.index,self.out_dir,self.fq_length,self.sd,self.thread,fq_cmd)
            python_tools.circ_call_process(cmd)
        quant_log = os.path.join(self.out_dir,'quant.cmd.log')
        python_tools.write_obj_to_file(cmd,quant_log)
        return cmd
    
    def kallisto_merge(self) :
        gene_quant_dict = {}
        for each_dir in self.sample_list :
            each_abs_path = os.path.join(self.out_dir,each_dir) 
            each_quant_file = os.path.join(each_abs_path,'abundance.tsv')
            if os.path.isdir(each_abs_path) :
                if not os.path.isfile(each_quant_file) :
                    sys.exit('sample :"%s" quant file is not exist !' % each_dir)
                quant_info = open(each_quant_file,'r') 
                for n,eachline in enumerate(quant_info) :
                    if n == 0 :
                        continue
                    else :
                        each_quant_info = eachline.strip().split('\t')
                        gene_id = each_quant_info[0]
                        est_counts = each_quant_info[3]
                        tpm = each_quant_info[4]
                        gene_quant_dict.setdefault(gene_id,{})[each_dir] = [est_counts,tpm]

        header = 'Transcript_ID\t%s\n' % '\t'.join(self.sample_list)
        est_counts_file = open(os.path.join(self.out_dir,'est_counts.xls'),'w')
        tpm_file = open(os.path.join(self.out_dir,'tpm.xls'),'w')
        est_counts_file.write(header)
        tpm_file.write(header)
        for each_gene in gene_quant_dict :
            each_gene_count_list = []
            each_gene_tpm_list   = []
            for each_sample in self.sample_list :
                each_gene_count_list.append(gene_quant_dict[each_gene][each_sample][0])
                each_gene_tpm_list.append(gene_quant_dict[each_gene][each_sample][1])
            est_counts_file.write('%s\t%s\n' % (each_gene,'\t'.join(each_gene_count_list)))
            tpm_file.write('%s\t%s\n' % (each_gene,'\t'.join(each_gene_tpm_list)))
        est_counts_file.close()
        tpm_file.close()        

class GO_enrich :
    def __init__(self) :
        self.transcript_dict = {}   
        self.gene_dict = {}
        self.gtf = None     
        self.target_length = None
        self.target_list = None
        self.go = None
        self.target_type = 'gene'
        self.out_dir = None
        self.auto_run = True
        self.plot_data = None
        self.compare_name = None
        self.enrich_dict = {}

    def get_target_length_table(self):
        if os.path.isfile(self.target_length) and os.stat(self.target_length).st_size :
            pass
        else :
            target_length_info = open(self.target_length,'w')
            if not self.transcript_dict :
                self.transcript_dict = get_transcript_info(self.gtf)
            if self.target_type == 'transcript' :
                for each_tr in self.transcript_dict :
                    each_tr_len = self.transcript_dict[each_tr]['length']
                    target_length_info.write('%s\t%s\n' % (each_tr,each_tr_len))
            else :
                self.gene_dict = get_gene_info(self.transcript_dict)
                for each_gene in self.gene_dict :
                    each_gene_len = python_tools.Median(self.gene_dict[each_gene]['transcript_len'])
                    target_length_info.write('%s\t%s\n' % (each_gene,each_gene_len))
            target_length_info.close()

    def add_enrich_info(self,compare_name,reg_flag,difffile) :
        self.enrich_dict.setdefault(compare_name,{})[reg_flag] = [difffile]

    def generate_goseq(self,difffile,output) :
        cmd = '%s/Rscript %s %s %s %s %s' % (R_BIN,GOSEQ,difffile,self.target_length,self.go,output)
        return cmd

    def run_GO_enrich(self) :
        cmd_list = []
        for each_compare in self.enrich_dict :
            for each_reg in self.enrich_dict[each_compare] :
                each_compare_out_dir = os.path.join(self.out_dir,each_compare)
                python_tools.circ_mkdir_unix(each_compare_out_dir)
                go_output = os.path.join(each_compare_out_dir,'%s.%s.GO.enrich.xls' % (each_compare,each_reg))
                self.enrich_dict[each_compare][each_reg].append(go_output)
                each_diff_file = self.enrich_dict[each_compare][each_reg][0]
                goseq_cmd = self.generate_goseq(each_diff_file,go_output)
                go_plot_flag = '%s.GO' % each_reg
                simple_bar_plot = enrich_bar_plot(each_compare,go_output,go_plot_flag,each_compare_out_dir)
                if self.auto_run :
                    python_tools.circ_call_process(goseq_cmd)
                    if os.path.isfile(go_output):
                        python_tools.circ_call_process(simple_bar_plot)
                cmd_list.extend([goseq_cmd, simple_bar_plot])
        return cmd_list

    def run_GO_DAG(self, gene_go_map):
        cmd_list = []
        for each_compare in self.enrich_dict:
            each_compare_out_dir = os.path.join(self.out_dir,each_compare)
            each_DAG_dir = os.path.join(each_compare_out_dir,'GO_DAG')
            cmd_list.append('mkdir -p %s' % each_DAG_dir)
            for each_reg in self.enrich_dict[each_compare]:
                go_output = os.path.join(each_compare_out_dir,'%s.%s.GO.enrich.xls' % (each_compare,each_reg))
                each_diff_file = self.enrich_dict[each_compare][each_reg][0]
                topgo_cmd = 'Rscript %s %s %s %s %s %s' % (TOPGO, gene_go_map, each_diff_file, go_output, each_reg, each_DAG_dir)
                if self.auto_run and os.path.exists(go_output):
                    python_tools.circ_mkdir_unix(each_DAG_dir)                                        
                    python_tools.circ_call_process(topgo_cmd)
                    cmd_list.append(topgo_cmd)
                else:
                    cmd_list.append(topgo_cmd)         
        return cmd_list            

class KEGG_enrich :
    def __init__(self) :
        self.seq = None
        self.all_blast_out = None
        self.species = None        
        self.ko_seq = None  
        self.auto_run = True
        self.blast_dict = {}
        self.target_type = 'gene'
        self.gtf = None
        self.transcript_dict = {} 
        self.enrich_dict = {}
        self.out_dir = None    
        self.sRNA_target = None     

    def check_blast_database(self) :
        suffix_list = ['phr','pin','psq']
        for each in suffix_list :
            if not os.path.isfile('%s.%s' % (self.ko_seq,each)) :
                os.system('%s/formatdb -i %s -p T' % (BLAST_BIN,self.ko_seq))
                break 

    def check_blast_program(self):
        seq_file_name = os.path.basename(self.seq)
        if 'pep' in seq_file_name :
            blast_program = 'blastp'
        elif 'cds' in seq_file_name :
            blast_program = 'blastx'
        else :
            sys.exit('Do not know use which blast program !')
        return blast_program   

    def get_blast_out(self):
        blast_bin = BLAST_BIN
        self.ko_seq = os.path.join(KO_PEP_DIR,'%s.pep.fasta' % self.species)
        self.check_blast_database()
        blast_pr = self.check_blast_program()
        cmd = '{blast_bin}/{blast_pr}  -query {self.seq} -db {self.ko_seq} -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -num_threads 1 -out {self.all_blast_out}'.format(**locals())
        if self.auto_run :
            python_tools.circ_call_process(cmd)
        return cmd

    def convert_blast(self):
        if self.target_type == 'gene' :
            if not self.transcript_dict :
                self.transcript_dict = get_transcript_info(self.gtf)

        with open(self.all_blast_out,'r') as kegg_blast_info :
            for eachline in kegg_blast_info :
                each_tr_blast = eachline.strip().split('\t')
                query_id = each_tr_blast[0]
                if self.target_type == 'gene' :
                    query_id = self.transcript_dict[each_tr_blast[0]]['gene_id']

                if query_id not in self.blast_dict :
                    self.blast_dict[query_id] = each_tr_blast[1:]
                else :
                    if self.blast_dict[query_id][-2] == each_tr_blast[-2] and self.blast_dict[query_id][-1] < each_tr_blast[-1]:
                        self.blast_dict[query_id] = each_tr_blast[1:]
                    elif self.blast_dict[query_id][-2] > each_tr_blast[-2] :
                        self.blast_dict[query_id] = each_tr_blast[1:]
                    else :
                        continue 

    def add_enrich_info(self,compare_name,reg_flag,diff_list) :
        self.enrich_dict.setdefault(compare_name,{})[reg_flag] = [diff_list]

    def check_KOBAS_out(self,kobas_out) :
        kobas_out_info = open(kobas_out,'r').readlines()
        flag = True
        for eachline in kobas_out_info :
            if not eachline.startswith("#") and len(eachline.strip().split('\t')) == 9 :
                flag = True
                break
        else :
            flag = False
        return flag

    def treat_KEGG_table(self,kegg_output):
        kegg_out_dir,kegg_out_name = os.path.split(kegg_output)
        kegg_tmp_file = os.path.join(kegg_out_dir,'tmp.%s' % kegg_out_name)
        os.system('mv %s %s' % (kegg_output,kegg_tmp_file))
        if self.check_KOBAS_out(kegg_tmp_file) :
            kegg_out_info = open(kegg_output,'w')
            with open(kegg_tmp_file,'r') as kegg_tmp_file_info :
                count = 0
                for eachline in kegg_tmp_file_info :            
                    if len(eachline.strip().split('\t')) == 9 :
                        if count == 0 and eachline.startswith("#") :
                            kegg_out_info.write(eachline)
                            count += 1
                        elif not eachline.startswith("#") :
                            kegg_out_info.write(eachline)
            kegg_out_info.close()
        os.system('rm %s' % (kegg_tmp_file))

    def generate_kobas(self,each_blast_out,kegg_output) :
        cmd = 'run_kobas.py -i {each_blast_out}  -t blastout:tab -s {self.species} -d K -o {kegg_output}'.format(**locals())
        return cmd

    def run_KEGG_enrich_old(self):
        cmd_list = []
        for each_compare in self.enrich_dict :
            for each_reg in self.enrich_dict[each_compare] :
                each_diff_list = self.enrich_dict[each_compare][each_reg][0]
                each_blast_out = os.path.join(self.out_dir,'%s.%s.blasttab' % (each_compare,each_reg))
                with open(each_blast_out,'w') as diff_blast_out :
                    for each_gene in each_diff_list :
                        if each_gene in self.blast_dict :
                            blast_info = '\t'.join(self.blast_dict[each_gene])
                            diff_blast_out.write('%s\t%s\n' % (each_gene,blast_info))

                each_compare_out_dir = os.path.join(self.out_dir,each_compare)
                python_tools.circ_mkdir_unix(each_compare_out_dir)
                kegg_output = os.path.join(each_compare_out_dir,'%s.%s.KEGG.enrich.xls' % (each_compare,each_reg))
                kegg_cmd = self.generate_kobas(each_blast_out,kegg_output)
                plot_name_tag = '%s.KEGG' % each_reg
                plot_cmd = enrich_bar_plot(each_compare,kegg_output,plot_name_tag,each_compare_out_dir)
                if self.auto_run :
                    python_tools.circ_call_process(kegg_cmd)
                    self.treat_KEGG_table(kegg_output)
                    cmd_list.append(kegg_cmd)
                    if os.path.isfile(kegg_output) :        
                        self.enrich_dict[each_compare][each_reg].append(kegg_output)                
                        self.enrich_dict[each_compare][each_reg].append(kegg_output)
                        python_tools.circ_call_process(plot_cmd)
                        cmd_list.append(plot_cmd)
                else:
                    treat_KEGG_table_cmd = 'python ~/scripts/treat_KEGG_table.py %s' % kegg_output
                    cmd_list.append(kegg_cmd)
                    cmd_list.append(treat_KEGG_table_cmd)
                    cmd_list.append(plot_cmd)
                    self.enrich_dict[each_compare][each_reg].append(kegg_output)                
        return cmd_list

    def run_KEGG_enrich(self):
        cmd_list = []
        blast_out_dir = os.path.join(self.out_dir, 'blast_out')
        python_tools.circ_mkdir_unix(blast_out_dir)
        for each_compare in self.enrich_dict :
            for each_reg in self.enrich_dict[each_compare] :
                each_diff_file = self.enrich_dict[each_compare][each_reg][0]
                each_compare_out_dir = os.path.join(self.out_dir,each_compare)
                python_tools.circ_mkdir_unix(each_compare_out_dir)
                kegg_output = os.path.join(each_compare_out_dir,'%s.%s.KEGG.enrich.xls' % (each_compare,each_reg))
                each_blast_out = os.path.join(blast_out_dir,'%s.%s.blasttab' % (each_compare,each_reg))
                extract_each_blast_cmd = 'python %s --id %s --table %s --output %s' % (EXTRACT_INF_BY_ID, each_diff_file, self.all_blast_out, each_blast_out)
                kegg_cmd = self.generate_kobas(each_blast_out, kegg_output)
                plot_name_tag = '%s.KEGG' % each_reg
                plot_cmd = enrich_bar_plot(each_compare,kegg_output,plot_name_tag,each_compare_out_dir)
                if self.auto_run :
                    python_tools.circ_call_process(extract_each_blast_cmd)
                    python_tools.circ_call_process(kegg_cmd)
                    self.treat_KEGG_table(kegg_output)
                    cmd_list.extend([extract_each_blast_cmd, kegg_cmd])
                    if os.path.isfile(kegg_output) :        
                        self.enrich_dict[each_compare][each_reg].append(kegg_output)                
                        self.enrich_dict[each_compare][each_reg].append(kegg_output)
                        python_tools.circ_call_process(plot_cmd)
                        cmd_list.append(plot_cmd)
                else:
                    treat_KEGG_table_cmd = 'python %s %s' % (TREAT_KEGG_TABLE, kegg_output)                    
                    cmd_list.extend([extract_each_blast_cmd, kegg_cmd, treat_KEGG_table_cmd, plot_cmd])
                    self.enrich_dict[each_compare][each_reg].append(kegg_output)                
        return cmd_list

    def run_kegg_pathview(self, diffdir, diffmethod, sRNA = False):
        cmd_list = []
        pathway_log_dir = os.path.join(self.out_dir,'kegg_pathway_logs')
        python_tools.circ_mkdir_unix(pathway_log_dir)
        for each_compare in self.enrich_dict:
            for each_reg in self.enrich_dict[each_compare]:
                each_compare_out_dir = os.path.join(self.out_dir,each_compare)
                # each_blast_out = os.path.join(self.out_dir,'%s.%s.blasttab' % (each_compare,each_reg))
                kegg_output = os.path.join(each_compare_out_dir,'%s.%s.KEGG.enrich.xls' % (each_compare,each_reg))
                diffout = glob.glob(r'%s/*%s.*.DE_results' % (diffdir,each_compare))[0]
                pathway_outdir = os.path.join(each_compare_out_dir,'%s_pathway' % each_reg)                
                pathview_check_log_file = os.path.join(pathway_log_dir,'%s.%s.log' % (each_compare,each_reg))
                if not sRNA:
                    pathview_cmd = 'python %s --kegg_table %s --blast_out %s --species %s --diff_out %s --diff_method %s --out_dir %s' % (PATHVIEW, kegg_output, self.all_blast_out, self.species, diffout, diffmethod, pathway_outdir)
                else:
                    pathview_cmd = 'python %s --kegg_table %s --blast_out %s --species %s --diff_out %s --diff_method %s --out_dir %s --gtf %s --sRNA_target %s' % (PATHVIEW_SRNA, kegg_output, self.all_blast_out, self.species, diffout, diffmethod, pathway_outdir, self.gtf, self.sRNA_target)                    
                pathview_check_cmd = 'python %s --kegg_table %s --pathway_dir %s --log_file %s' % (PATHVIEW_CK, kegg_output, pathway_outdir, pathview_check_log_file)
                if self.auto_run and os.path.exists(kegg_output): 
                    python_tools.circ_mkdir_unix(pathway_outdir)                                         
                    python_tools.circ_call_process(pathview_cmd)
                    python_tools.circ_call_process(pathview_check_cmd)
                    cmd_list.extend([pathview_cmd,pathview_check_cmd])
                else:
                    mkdir_cmd = 'mkdir -p %s' % pathway_outdir
                    cmd_list.extend([mkdir_cmd, pathview_cmd,pathview_check_cmd])
        return cmd_list

def get_ym_md5(md5_file):
    md5_file_info = open(md5_file).readlines()
    md5 = md5_file_info[0].split()[0]
    return md5


class rawdata :
    def __init__(self) :
        self.name = None
        self.read1 = []
        self.read2 = []
        self.md5 = []
        self.pre_num = 0

    def check_rawdata(self):
        assert len(self.read1) == len(self.read2)

    def cp_rawdata_to_custum(self,out_dir) :
        cmd = '\necho "#### cp rawdata of sample : %s start ####"' % self.name
        cmd = [cmd]      
        each_sample_dir = os.path.join(out_dir,self.name)
        if not os.path.exists(each_sample_dir) :
            mkdir_cmd = 'mkdir "%s"' % each_sample_dir
            cmd.append(mkdir_cmd)

        self.check_rawdata()
        if len(self.read1) == 0 :
            sys.exit('%s without rawdata !' % self.name)
        elif len(self.read1) == 1 :            
            r1_output = os.path.join(each_sample_dir,'%s_R1.fastq.gz' % (self.name))
            r2_output = os.path.join(each_sample_dir,'%s_R2.fastq.gz' % (self.name))
            cmd1 = 'cp "%s" "%s"' % (self.read1[0],r1_output)
            cmd2 = 'cp "%s" "%s"' % (self.read2[0],r2_output)
            if os.path.isfile(r1_output) and os.stat(self.read1[0]).st_size == os.stat(r1_output).st_size :   
                cmd.append('echo "#### %s already exists! ####"' % r1_output)
            else :
                cmd.append(cmd1)
            if os.path.isfile(r2_output) and os.stat(self.read2[0]).st_size == os.stat(r2_output).st_size :
                cmd.append('echo "#### %s already exists! ####"' % r2_output)
            else :
                cmd.append(cmd2)
        else :  
            for n in xrange(len(self.read1)) :
                r1_output = os.path.join(each_sample_dir,'%s_%s_R1.fastq.gz' % (self.name,n+1))
                r2_output = os.path.join(each_sample_dir,'%s_%s_R2.fastq.gz' % (self.name,n+1))
                cmd1 = 'cp "%s" "%s"' % (self.read1[n],r1_output)
                cmd2 = 'cp "%s" "%s"' % (self.read2[n],r2_output)
                if os.path.isfile(r1_output) and os.stat(self.read1[n]).st_size == os.stat(r1_output).st_size :   
                    cmd.append('echo "#### %s already exists! ####"' % r1_output)
                else:
                    cmd.append(cmd1)
                if os.path.isfile(r2_output) and os.stat(self.read2[n]).st_size == os.stat(r2_output).st_size :
                    cmd.append('echo "#### %s already exists! ####"' % r2_output)
                else:
                    cmd.append(cmd2)
        cmd_end = 'echo "#### cp rawdata of sample : %s end ####"' % self.name
        cmd.append(cmd_end)
        return cmd

    def merge_rawdata(self,out_dir) :
        cmd = '\necho "#### get rawdata of sample : %s start ####"' % self.name
        cmd = [cmd]  
        self.check_rawdata()
        r1_output = os.path.join(out_dir,'%s_R1.fastq.gz' % (self.name))
        r2_output = os.path.join(out_dir,'%s_R2.fastq.gz' % (self.name))
        if len(self.read1) == 0 :
            sys.exit('%s without rawdata !' % self.name)
        elif len(self.read1) == 1 :            
            cmd1 = 'cp "%s" "%s"' % (self.read1[0],r1_output)
            cmd2 = 'cp "%s" "%s"' % (self.read2[0],r2_output)
            if os.path.isfile(r1_output) and os.stat(self.read1[0]).st_size == os.stat(r1_output).st_size :   
                cmd.append('echo "#### %s already exists! ####"' % r1_output)
            else :
                cmd.append(cmd1)
            if os.path.isfile(r2_output) and os.stat(self.read2[0]).st_size == os.stat(r2_output).st_size :
                cmd.append('echo "#### %s already exists! ####"' % r2_output)  
            else :
                cmd.append(cmd2)                    
        else :
            r1_files = ' '.join(['"%s"' % each for each in self.read1])
            r2_files = ' '.join(['"%s"' % each for each in self.read2])
            cmd1 = 'zcat %s | gzip > "%s"' % (r1_files,r1_output)
            cmd2 = 'zcat %s | gzip > "%s"' % (r2_files,r2_output)
            if len(self.read1) > self.pre_num:
                cmd.extend([cmd1,cmd2])
            elif len(self.read1) == self.pre_num:
                if os.path.isfile(r1_output) and os.stat(r1_output).st_size :
                    cmd.append('echo "#### %s already exists! ####"' % r1_output)
                else :
                    cmd.append(cmd1)
                if os.path.isfile(r2_output) and os.stat(r2_output).st_size :
                    cmd.append('echo "#### %s already exists! ####"' % r2_output)  
                else :
                    cmd.append(cmd2)                 
            else:
                sys.exit('less fastq file than last time!')
        cmd_end = 'echo "#### get rawdata of sample : %s end ####"' % self.name
        cmd.append(cmd_end)
        return cmd

    def add_sup_data(self,out_dir) :
        cmd = '\necho "#### cp rawdata of sample : %s start ####"' % self.name
        cmd = [cmd]      
        each_sample_dir = os.path.join(out_dir,self.name)
        if not os.path.exists(each_sample_dir) :
            mkdir_cmd = 'mkdir %s' % each_sample_dir
            cmd.append(mkdir_cmd)

        self.check_rawdata()
        if len(self.read1) == 0 :
            sys.exit('%s without rawdata !' % self.name)
        elif len(self.read1) == 1 :            
            r1_output = os.path.join(each_sample_dir,'%s_2_R1.fastq.gz' % (self.name))
            r2_output = os.path.join(each_sample_dir,'%s_2_R2.fastq.gz' % (self.name))
            cmd1 = 'cp "%s" "%s"' % (self.read1[0],r1_output)
            cmd2 = 'cp "%s" "%s"' % (self.read2[0],r2_output)
            if os.path.isfile(r1_output) and os.stat(self.read1[0]).st_size == os.stat(r1_output).st_size :   
                cmd.append('echo "#### %s already exists! ####"' % r1_output)
            else :
                cmd.append(cmd1)
            if os.path.isfile(r2_output) and os.stat(self.read2[0]).st_size == os.stat(r2_output).st_size :
                cmd.append('echo "#### %s already exists! ####"' % r2_output)
            else :
                cmd.append(cmd2)
        else :  
            for n in xrange(len(self.read1)) :
                r1_output = os.path.join(each_sample_dir,'%s_%s_R1.fastq.gz' % (self.name,n+2))
                r2_output = os.path.join(each_sample_dir,'%s_%s_R2.fastq.gz' % (self.name,n+2))
                cmd1 = 'cp "%s" "%s"' % (self.read1[n],r1_output)
                cmd2 = 'cp "%s" "%s"' % (self.read2[n],r2_output)
                if os.path.isfile(r1_output) and os.stat(self.read1[n]).st_size == os.stat(r1_output).st_size :   
                    cmd.append('echo "#### %s already exists! ####"' % r1_output)
                else:
                    cmd.append(cmd1)
                if os.path.isfile(r2_output) and os.stat(self.read2[n]).st_size == os.stat(r2_output).st_size :
                    cmd.append('echo "#### %s already exists! ####"' % r2_output)
                else:
                    cmd.append(cmd2)
        cmd_end = 'echo "#### cp rawdata of sample : %s end ####"' % self.name
        cmd.append(cmd_end)
        return cmd

    def check_dna_md5(self, out_dir):
        log = []
        each_sample_dir = os.path.join(out_dir,self.name)
        self.check_rawdata()
        if len(self.read1) == 0:
            sys.exit('%s without rawdata !' % self.name)
        else:
            for n in xrange(len(self.read1)):
                for m,each_read_file in enumerate([self.read1[n],self.read2[n]]):
                    read_num = m+1
                    seq_num = n+1                    
                    each_seq_md5_file = '%s.md5' % each_read_file
                    if os.path.isfile(each_seq_md5_file):
                        each_seq_md5 = get_ym_md5(each_seq_md5_file)
                    else:
                        log.append('%s not exists!' % each_seq_md5_file)
                        continue                    
                    if len(self.read1) == 1:
                        each_read_name = '%s_R%s.fastq.gz' % (self.name,read_num)
                        each_read_output = os.path.join(each_sample_dir, each_read_name)
                    else:
                        each_read_name = '%s_%s_R%s.fastq.gz' % (self.name,seq_num,read_num)
                        each_read_output = os.path.join(each_sample_dir, each_read_name)
                    if os.path.isfile(each_read_output):
                        each_custumer_md5 = md5sum(each_read_output)
                        if each_custumer_md5 == each_seq_md5:
                            log.append('%s md5 check ok!' % each_read_output)
                            self.md5.append('%s\t%s' % (each_read_name,each_custumer_md5))
                        else:
                            each_disk_md5 = md5sum(each_read_file)
                            if each_disk_md5 == each_seq_md5:
                                log.append('%s cp error!' % each_read_output)
                            else:
                                log.append('%s md5 value not consists with %s!' % (each_read_file, each_seq_md5_file))
                    else:
                        log.append('%s not exists!' % each_read_output)
        return log

    def check_rna_md5(self, out_dir):    
        log = []
        if len(self.read1) == 0:
            log.append('%s without rawdata !' % self.name)
        else:
            for m,each_read_file in enumerate([self.read1[-1],self.read2[-1]]):
                read_num = m+1
                each_read_name = '%s_R%s.fastq.gz' % (self.name,read_num)
                each_read_output = os.path.join(out_dir,each_read_name)
                each_md5_output = os.path.join(out_dir,'%s.md5' % (each_read_name))                             
                if len(self.read1) == 1:
                    each_seq_md5_file = '%s.md5' % each_read_file
                    if os.path.isfile(each_seq_md5_file):
                        each_seq_md5 = get_ym_md5(each_seq_md5_file)
                    else:
                        log.append('%s not exists!' % each_seq_md5_file)
                        continue
                    if os.path.exists(each_read_output):
                        each_read_output_md5 = md5sum(each_read_output)
                    else:
                        log.append('%s not exists!' % each_read_output)
                        continue
                    if each_read_output_md5 == each_seq_md5:
                        log.append('%s md5 check ok!' % each_read_output)
                        self.md5.append('%s\t%s' % (each_read_name,each_read_output_md5))
                    else:
                        each_read_md5 = md5sum(each_read_file)
                        if each_read_md5 == each_seq_md5:
                            log.append('%s cp error!' % each_read_output)
                        else:
                            log.append('%s md5 value not consists with %s!' % (each_read_file, each_seq_md5_file))
                else:
                    each_read_output_last_line = os.popen('zcat "%s" | tail -4' % each_read_output).read()
                    each_read_file_last_line = os.popen('zcat "%s" | tail -4' % each_read_file).read()
                    if each_read_output_last_line == each_read_file_last_line:
                        each_read_output_md5 = md5sum(each_read_output)
                        self.md5.append('%s\t%s' % (each_read_name,each_read_output_md5))
                    else:
                        log.append('%s data zcat error!' % self.name)
        return log    

    def get_dna_md5(self):
        log = []
        self.check_rawdata()
        if len(self.read1) == 0:
            sys.exit('%s without rawdata !' % self.name)
        else:
            for n in xrange(len(self.read1)):
                for m,each_read_file in enumerate([self.read1[n],self.read2[n]]):
                    read_num = m+1
                    seq_num = n+1                    
                    each_seq_md5_file = '%s.md5' % each_read_file
                    each_seq_size = round(float(os.stat(each_read_file).st_size)/1024**3, 2)
                    if os.path.isfile(each_seq_md5_file):
                        each_seq_md5 = get_ym_md5(each_seq_md5_file)
                    else:
                        log.append('%s not exists!' % each_seq_md5_file)
                        continue                    
                    if len(self.read1) == 1:
                        each_read_name = '%s_R%s.fastq.gz' % (self.name,read_num)
                    else:
                        each_read_name = '%s_%s_R%s.fastq.gz' % (self.name,seq_num,read_num)
                    self.md5.append('%s\t%s\t%sG\t%s' % (self.name, each_read_name, each_seq_size, each_seq_md5))
        return log

class make_maf :
    def __init__(self) :
        self.target_species = None
        self.query_species = []
        self.target_species_chrom = []
        self.subfa_file = None
        self.subfa_dict = {}
        self.lastzNear="C=0 E=150 H=0 K=4500 L=3000 M=254 O=600 Q=/public/database/ucsc/conservation/human_chimp.v2.q T=2 Y=15000"
        self.lastzMedium="C=0 E=30 H=0 K=3000 L=3000 M=50 O=400 T=1 Y=9400"
        self.lastzFar="C=0 E=30 H=2000 K=2200 L=6000 M=50 O=400 Q=/public/database/ucsc/conservation/HoxD55.q T=2 Y=3400"
        self.out_dir = None
        self.tree = None

    def get_subfa_dict(self) :
        with open(self.subfa_file) as subfa_file_info :
            for eachline in subfa_file_info :
                eachline_info = eachline.strip().split('\t')
                species_id = eachline_info[0]
                subfa_dir = eachline_info[1]
                if species_id in self.query_species or species_id == self.target_species :
                    self.subfa_dict[species_id] = subfa_dir

    def pairwise_alignment(self):
        run_script_list = []
        for each_sp in self.query_species :
            each_cmp_name = '%s_vs_%s' % (self.target_species,each_sp)
            each_sp_outdir = os.path.join(self.out_dir, each_cmp_name )
            each_sp_script_dir = os.path.join(each_sp_outdir,'scripts')
            each_sp_result_dir = os.path.join(each_sp_outdir,'results')
            each_sp_result_split_dir = os.path.join(each_sp_result_dir,'split')
            map(python_tools.circ_mkdir_unix,[each_sp_script_dir,each_sp_result_split_dir])
            each_sp_subfa = self.subfa_dict[each_sp]
            target_species_subfa = self.subfa_dict[self.target_species]
            each_sp_subfa_list = os.listdir(each_sp_subfa)
            each_target_species_subfa = os.listdir(target_species_subfa)
            combined_aligment_scirpt = os.path.join(each_sp_script_dir,'%s.aligment.sh' % each_cmp_name)
            combined_aligment_scirpt_list = []
            combined_merge_script = os.path.join(each_sp_script_dir,'%s.merge.sh' % each_cmp_name)
            combined_merge_script_list = []  
            python_tools.write_obj_to_file('#! /bin/bash',combined_merge_script)          
            for each_target in each_target_species_subfa :
                for each_query in each_sp_subfa_list :
                    each_target_file = os.path.join(target_species_subfa,each_target)
                    each_query_file = os.path.join(each_sp_subfa,each_query)
                    each_cmd = os.path.join(each_sp_script_dir,'%s_vs_%s.sh' % (each_target,each_query))
                    each_out = os.path.join(each_sp_result_split_dir,'%s_vs_%s.maf' % (each_target,each_query))
                    each_cmd_info =  'echo "####%s align %s begins####"\n' % (each_target,each_query)
                    each_cmd_info += 'date\n'
                    #each_cmd_info.write('lastz {each_target_file} {each_query_file} --notransition --step=20 --ambiguous=iupac {self.lastzMedium} --format=maf > {each_out}\n'.format(**locals()))
                    each_cmd_info += 'lastz {each_target_file} {each_query_file} --notransition --step=20 --ambiguous=iupac {self.lastzMedium} --format=maf > {each_out}\n'.format(**locals())
                    each_cmd_info += 'date\n'
                    each_cmd_info +=  'echo "####%s align %s ends####"\n' % (each_target,each_query)
                    python_tools.write_obj_to_file(each_cmd_info,each_cmd)
                    combined_aligment_scirpt_list.append(each_cmd)
                merge_info = 'cat %s/%s* > %s/%s_vs_%s.maf' % (each_sp_result_split_dir,each_target,each_sp_result_dir,each_target,each_sp)
                combined_merge_script_list.append(merge_info)
            python_tools.multi_process_shell_script(combined_aligment_scirpt_list,combined_aligment_scirpt,12)
            python_tools.write_obj_to_file(combined_merge_script_list,combined_merge_script,True)
            merge_cmd_line = '\nwait\necho "#### merge aligment results ####"\nsh %s\necho "#### merge finished ####"\n' % combined_merge_script
            python_tools.write_obj_to_file(merge_cmd_line,combined_aligment_scirpt,True)
            os.system('chmod +x %s' % combined_aligment_scirpt)
            os.system('chmod +x %s' % combined_merge_script)
            run_script_list.append(combined_aligment_scirpt)
            # python_tools.circ_call_process(combined_aligment_scirpt)
        return run_script_list

    def merge_maf(self):
        merge_script_list = []
        merge_dir = os.path.join(self.out_dir,'merge_maf')
        for each_chrom in self.target_species_chrom:
            each_chrom_dir = os.path.join(merge_dir,each_chrom)
            each_tmp_dir = os.path.join(each_chrom_dir,'tmp')
            merge_script_list.append('echo "#### start merge chromosome %s ####"' % each_chrom)
            if not os.path.exists(each_chrom_dir):
                merge_script_list.append('mkdir -p %s' % each_tmp_dir)
            merge_script_list.append('cd %s' % each_chrom_dir)
            for each_sp in self.query_species:
                each_cmp_name = '%s_vs_%s' % (self.target_species,each_sp)
                each_sp_outdir = os.path.join(self.out_dir, each_cmp_name )
                each_sp_result_dir = os.path.join(each_sp_outdir,'results')
                each_out = os.path.join(each_sp_result_dir,'%s.%s.fa_vs_%s.maf' % (self.target_species,each_chrom,each_sp))
                each_merge_name = os.path.join(each_chrom_dir,'%s.%s.sing.maf' % (self.target_species,each_sp))
                link_cmd = 'ln -s %s %s' % (each_out,each_merge_name)
                merge_script_list.append(link_cmd)
            merge_script_list.append('roast E=%s T=%s "%s" %s/*.*.sing.maf %s/%s.maf' % (self.target_species,each_tmp_dir,self.tree,each_chrom_dir,merge_dir,each_chrom))
            merge_script_list.append('rm -r %s' % each_tmp_dir)
            merge_script_list.append('echo "#### finish merge chromosome %s ####"' % each_chrom)
        return merge_script_list        

def get_sample_name(fq_data):
    match = re.search(r'(\S+)_([1,2]).clean.fq.gz',fq_data)
    if match:
        fq_data_name = match.groups()[0]
        read_num = match.groups()[1]
        return fq_data_name,read_num
    else:
        sys.exit('Wrong fq file name! : %s' % fq_data)

def launch_jobs(jobfiles,out_dir,platform = 'th',name = "job"):
    launch_job_file = os.path.join(out_dir,'launch_%s.sh' % name)
    launch_job_cmds = []
    if isinstance(jobfiles,str):
        jobfiles = [jobfiles]
    for each_job_file in jobfiles:
        if os.path.isfile(each_job_file):
            if platform == 'th':
                os.system('chmod +x %s' % each_job_file)
                cmd_line = 'yhrun -n 1 -N 1 -c 24 -x cn[9531] -p work %s &' % each_job_file                
            else:
                cmd_line = 'sh %s\nwait' % each_job_file
            launch_job_cmds.append(cmd_line)
        else:
            sys.exit('%s not exits!' % each_job_file)
    python_tools.write_obj_to_file(launch_job_cmds,launch_job_file)


def gtf_length_filter(gtf, tr_dict, length = 200):
    out_list = []
    for eachline in GFF_Reader(gtf):
        tr_id = eachline.attr['transcript_id']
        try:
            if tr_dict[tr_id]['length'] > 200:
                out_list.append('%s;' % eachline.get_gff_line().strip())
        except:
            print tr_id
            sys.exit(0)
    return out_list

class RNAseq_pipeline():
    def __init__(self):
        ## pipeline type
        # self.pipe = 'mRNA'
        # self.mRNA_analysis_list = ['quant', 'diff', 'enrich', ]
        # self.lncRNA_analysis_list = ['mapping', 'assembly', 'filter', 'quant', 'diff', 'target', 'enrich']              
        # self.pipe_content_dict = {
        #     'mRNA':self.mRNA_analysis_list, 
        #     'lncRNA':self.lncRNA_analysis_list,
        # }

        ## general parameters
        self.cleandata_dir = None
        self.sample_dict = {}
        ## qc
        self.qc_dir = ""
        self.qc_thread = 2
        ## mapping parameters
        self.mapping_thread = 8
        self.genome_bowtie2_index = None
        self.trans_bowtie2_index = None
        self.gtf = None
        self.mapping_program = 'tophat2'
        self.libtype = 'fr-unstranded'
        self.mapping_dir = None
        self.platform = 'th'
        ## assembly parameters
        self.assembly_thread = 8
        self.assembly_dir = None
        ## quant parameters
        self.transcript_fa = None
        self.gene_trans = None
        self.group_sample = None
        self.quant_program = 'kallisto'
        self.quant_thread = 4
        self.quant_dir = None
        self.qvalue = 0.05
        self.logfc = 1
        self.diff_program = 'edgeR'
        self.diff_dir = None
        ## enrichment analysis
        self.enrich_thread = 1
        self.compare_list = []
        self.kegg_species = None
        self.go = None
        self.topgo = None
        self.blast_out = None
        self.gene_length = None
        self.enrich_dir = None
        self.run_go = True
        self.run_kegg = True
        ## results generate
        self.anno_files = ""
        self.result_dir = None

    def get_cleandata(self):
        sample_dict_json_file = os.path.join(self.cleandata_dir,'sample_info.json')
        if os.path.isfile(sample_dict_json_file):
            self.sample_dict = python_tools.load_fn_to_obj(sample_dict_json_file)
        else:
            cleandata_list = os.listdir(self.cleandata_dir)
            for each_fq in cleandata_list:
                if each_fq.endswith('clean.fq.gz'):
                    fq_data_name,read_num = get_sample_name(each_fq)
                    each_fq_path = os.path.join(self.cleandata_dir,each_fq)
                    self.sample_dict.setdefault(fq_data_name,{})[read_num] = each_fq_path
            python_tools.write_obj_to_json(self.sample_dict, sample_dict_json_file)                

    def cmd_to_scripts(self, cmd_list, thread, out_dir, server = 'server228', name = 'job', split = True):
        parallel_num = SERVER_CORE_DICT[server]//thread
        script_num = 0
        script_cmd_list = ['#!/bin/sh\ndate\necho "#### start ####"']
        script_file_list = []
        if parallel_num >= len(cmd_list) or (not split):
            script_path = os.path.join(out_dir,'run_%s.sh' % name)
            script_cmd_list.extend(cmd_list) 
            script_cmd_list.append('wait\necho "#### finished ####"\ndate')
            python_tools.write_obj_to_file(script_cmd_list,script_path)    
            script_file_list = [script_path]
        else:
            for n,each_cmd in enumerate(cmd_list):
                script_cmd_list.append(each_cmd)                       
                if (n+1) % parallel_num == 0 or (n+1) == len(cmd_list):
                    script_num += 1
                    script_path = os.path.join(out_dir,'run_%s_part%s.sh' % (name, script_num))
                    script_file_list.append(script_path)
                    script_cmd_list.append('wait\necho "#### finished ####"\ndate')
                    python_tools.write_obj_to_file(script_cmd_list,script_path)                
                    script_cmd_list = ['#!/bin/sh\ndate\necho "#### start ####"']
        launch_jobs(script_file_list,out_dir,self.platform,name)

    def run_qc(self):
        if not self.cleandata_dir:
            sys.exit('cleandata dierectory not exists!')
        else:
            if not self.sample_dict:
                self.get_cleandata()
        qc_cmd_list = []
        python_tools.circ_mkdir_unix(self.qc_dir)
        for each_sample in self.sample_dict:
            read1 = self.sample_dict[each_sample]['1']
            read2 = self.sample_dict[each_sample]['2']            
            if self.mapping_program == 'tophat2':
                qc_cmd_list.append('fastqc {read1} {read2} -o {self.qc_dir}'.format(**locals()))
        self.cmd_to_scripts(qc_cmd_list, self.qc_thread, self.qc_dir, self.platform)


    def check_quant_index(self):
        if self.quant_program == 'kallisto':
            quant_index = '%s.kallisto_idx' % self.transcript_fa
            if os.path.isfile(quant_index):
                return True
        elif self.quant_program == 'salmon':
            quant_index = '%s.salmon_quasi.idx' % self.transcript_fa
            if os.path.isdir(quant_index):
                return True
        else:
            sys.exit('Not support quant program: %s' % self.quant_program)

    def tophat2_mapping(self, read1, read2, out_dir):
        tophat2_mapping_cmd = 'tophat --library-type {self.libtype} -p {self.mapping_thread} -G {self.gtf} --transcriptome-index={self.trans_bowtie2_index} -o {out_dir} {self.genome_bowtie2_index} {read1} {read2} &'.format(**locals())
        return tophat2_mapping_cmd

    def stringtie_assembly(self, bam, name, out_dir):
        stringtie_cmd = 'stringtie {bam} -p {self.assembly_thread} -G {self.gtf} -o {out_dir}/{name}.gtf &'.format(**locals())
        return stringtie_cmd

    def run_mapping(self):
        if not self.cleandata_dir:
            sys.exit('cleandata dierectory not exists!')
        else:
            if not self.sample_dict:
                self.get_cleandata()
        mapping_cmd_list = []
        for each_sample in self.sample_dict:
            read1 = self.sample_dict[each_sample]['1']
            read2 = self.sample_dict[each_sample]['2']
            each_sample_dir = os.path.join(self.mapping_dir,each_sample)
            python_tools.circ_mkdir_unix(each_sample_dir)
            if self.mapping_program == 'tophat2':
                mapping_cmd_list.append(self.tophat2_mapping(read1,read2,each_sample_dir))
        self.cmd_to_scripts(mapping_cmd_list,self.mapping_thread,self.mapping_dir, self.platform)

    def run_assembly(self):
        if not self.cleandata_dir:
            sys.exit('cleandata dierectory not exists!')
        else:
            if not self.sample_dict:
                self.get_cleandata()
        if not os.path.exists(self.assembly_dir):
            python_tools.circ_mkdir_unix(self.assembly_dir)
        assembly_cmd_list = []
        for each_sample in self.sample_dict:
            each_sample_dir = os.path.join(self.mapping_dir, each_sample)
            each_sample_bam = os.path.join(each_sample_dir, 'accepted_hits.bam')
            if os.path.isfile(each_sample_bam):
                assembly_cmd_list.append(self.stringtie_assembly(each_sample_bam, each_sample, self.assembly_dir))
            else:
                sys.exit('%s not exists!' % each_sample_bam)
        self.cmd_to_scripts(assembly_cmd_list, self.assembly_thread, self.assembly_dir, self.platform)

    def run_quant(self):
        if self.quant_program == 'salmon':
            os.system('export LC_ALL="C"')
        quant_script_dir = os.path.join(self.quant_dir,'scripts')
        python_tools.circ_mkdir_unix(quant_script_dir)
        if not self.cleandata_dir:
            sys.exit('cleandata dierectory not exists!')
        else:
            if not self.sample_dict:
                self.get_cleandata()
        # quant_cmd_list = ['date\necho "#### quanticification for each sample ####"']
        quant_cmd_list = []
        quant_results = []
        diff_cmd_list = []
        for n,each_sample in enumerate(self.sample_dict):
            read1 = self.sample_dict[each_sample]['1']
            read2 = self.sample_dict[each_sample]['2']
            each_sample_dir = os.path.join(self.quant_dir,each_sample)
            if self.quant_program == 'kallisto':
                each_quant_results = os.path.join(each_sample_dir,'abundance.tsv.genes')
            elif self.quant_program == 'salmon':
                each_quant_results = os.path.join(each_sample_dir,'quant.sf.genes')
            quant_results.append(each_quant_results)
            if n == 0 and not (self.check_quant_index()):
                index_script = ['#!/bin/sh\ndate\necho "#### build index start ####"']
                if self.libtype == 'fr-unstranded':
                    quant_cmd = 'align_and_estimate_abundance.pl --transcripts {self.transcript_fa} --prep_reference --seqType fq --left {read1} --right {read2} --est_method {self.quant_program} --gene_trans_map {self.gene_trans} --output_dir {each_sample_dir} \nwait'.format(**locals())
                elif self.libtype == 'fr-firststrand' and self.quant_program != 'kallisto':
                    quant_cmd = 'align_and_estimate_abundance.pl --transcripts {self.transcript_fa} --SS_lib_type RF --prep_reference --seqType fq --left {read1} --right {read2} --est_method {self.quant_program} --gene_trans_map {self.gene_trans} --output_dir {each_sample_dir} \nwait'.format(**locals())
                else:
                    sys.exit('wrong parameters combination! : libtype %s quant_program %s' % (self.libtype, self.quant_program))
                index_script.append(quant_cmd)
                index_script.append('date\necho "#### build index finished ####"')
                self.cmd_to_scripts(index_script,self.quant_thread, quant_script_dir, self.platform, 'index_build', False)    
            else:
                if self.libtype == 'fr-unstranded':
                    quant_cmd = 'align_and_estimate_abundance.pl --transcripts {self.transcript_fa} --seqType fq --left {read1} --right {read2} --est_method {self.quant_program} --gene_trans_map {self.gene_trans} --output_dir {each_sample_dir} &'.format(**locals())
                elif self.libtype == 'fr-firststrand' and self.quant_program != 'kallisto':
                    quant_cmd = 'align_and_estimate_abundance.pl --transcripts {self.transcript_fa} --SS_lib_type RF --seqType fq --left {read1} --right {read2} --est_method {self.quant_program} --gene_trans_map {self.gene_trans} --output_dir {each_sample_dir} &'.format(**locals())
                else:
                    sys.exit('wrong parameters combination! : libtype %s quant_program %s' % (self.libtype, self.quant_program))
                quant_cmd_list.append(quant_cmd)
        self.cmd_to_scripts(quant_cmd_list, self.quant_thread, quant_script_dir, self.platform, 'quant')
        diff_cmd_list.append('date\necho "#### Merge individual quant table ####"')
        quant_results_line = ' '.join(quant_results)
        diff_cmd_list.append('abundance_estimates_to_matrix.pl --est_method {self.quant_program} --out_prefix {self.quant_dir}/genes --name_sample_by_basedir {quant_results_line}'.format(**locals()))
        ## sample correlation
        diff_cmd_list.append('wait\ndate\n\necho "#### QC Samples and Replicates ####"')
        sample_qc_dir = os.path.join(self.quant_dir,'sample_correlation')
        python_tools.circ_mkdir_unix(sample_qc_dir)
        diff_cmd_list.append('cd %s' % sample_qc_dir)
        diff_cmd_list.append('PtR --matrix {self.quant_dir}/genes.counts.matrix --samples {self.group_sample} --log2 --compare_replicates --sample_cor_matrix'.format(**locals()))
        ## diff analysis
        python_tools.circ_mkdir_unix(self.diff_dir)
        diff_cmd_list.append('wait\ndate\n\necho "#### Differential analysis ####"')
        diff_cmd_list.append('cd %s' % self.diff_dir)
        diff_cmd_list.append('run_DE_analysis.pl --matrix {self.quant_dir}/genes.counts.matrix --method {self.diff_program} --samples_file {self.group_sample} --output {self.diff_dir}'.format(**locals()))
        diff_cmd_list.append('analyze_diff_expr.pl -P {self.qvalue} -C {self.logfc} --matrix {self.quant_dir}/genes.TMM.EXPR.matrix --samples {self.group_sample} --max_genes_clust 20000'.format(**locals()))
        for ptree in [10,15,20]:
            diff_cmd_list.append('define_clusters_by_cutting_tree.pl --Ptree {ptree} -R {self.diff_dir}/diffExpr.P{self.qvalue}_C{self.logfc}.matrix.RData'.format(**locals()))
        diff_cmd_list.append('echo "#### quanticification analysis finished ####"')
        self.cmd_to_scripts(diff_cmd_list, self.quant_thread, self.diff_dir, self.platform, 'diff', False)

    def get_compare_list(self):
        for each_file in os.listdir(self.diff_dir):
            if each_file.endswith('DE_results'):
                each_compare = re.search(r'genes.counts.matrix.(.*?).%s.DE_results' % self.diff_program,each_file).groups()[0]
                if each_compare not in self.compare_list:
                    self.compare_list.append(each_compare)

    def get_diff_gene(self, difffile, outgene):
        diff_gene_list = []
        if os.path.isfile(difffile):
            with open(difffile) as difffile_info:
                for n,eachline in enumerate(difffile_info):
                    if n!= 0:
                        eachline_info = eachline.split('\t')
                        diff_gene_list.append(eachline_info[0])
        else:
            sys.exit('%s not exists!' % difffile)
        python_tools.write_obj_to_file(diff_gene_list,outgene)

    def run_enrich(self):
        if not self.compare_list:
            self.get_compare_list()
        ## get GO parameters
        if self.run_go:
            my_go_enrich = GO_enrich()
            my_go_enrich.target_length = self.gene_length
            my_go_enrich.go = self.go
            my_go_enrich.out_dir = os.path.join(self.enrich_dir,'GO')   
            my_go_enrich.auto_run = False     
        ## get KEGG parameters
        if self.run_kegg:
            my_kegg_enrich = KEGG_enrich()
            my_kegg_enrich.out_dir = os.path.join(self.enrich_dir,'KEGG')
            diff_gene_list_dir = os.path.join(self.diff_dir,'Diff_list')
            my_kegg_enrich.all_blast_out = self.blast_out
            my_kegg_enrich.species = self.kegg_species
            my_kegg_enrich.auto_run = False
        for each_compare in self.compare_list:
            each_GO_compare_outdir = os.path.join(my_go_enrich.out_dir,each_compare)
            each_KEGG_compare_outdir = os.path.join(my_kegg_enrich.out_dir,each_compare)
            map(python_tools.circ_mkdir_unix,[each_GO_compare_outdir,each_KEGG_compare_outdir, diff_gene_list_dir])
            cond_list = each_compare.split('_vs_')
            diff_gene_file_list = []
            for each_cond in cond_list:
                each_name = '%s-UP' % each_cond
                each_cond_up_file = glob.glob(r'%s/*%s.*.DE_results.*.%s-UP.subset' % (self.diff_dir,each_compare,each_cond))[0]
                each_cond_up_gene = os.path.join(diff_gene_list_dir,'%s.%s-UP.list' % (each_compare, each_cond))
                self.get_diff_gene(each_cond_up_file, each_cond_up_gene)
                list_obj = [each.strip() for each in open(each_cond_up_gene)]
                if list_obj:
                    my_kegg_enrich.add_enrich_info(each_compare,each_name,each_cond_up_gene)  
                    my_go_enrich.add_enrich_info(each_compare,each_name,each_cond_up_gene)
                    diff_gene_file_list.append(each_cond_up_gene)
            all_diff_gene = os.path.join(diff_gene_list_dir,'%s.ALL.list' % each_compare)
            python_tools.merge_files(diff_gene_file_list,all_diff_gene)
            list_obj = [each.strip() for each in open(all_diff_gene)]
            if list_obj:
                my_kegg_enrich.add_enrich_info(each_compare,'ALL',all_diff_gene)  
                my_go_enrich.add_enrich_info(each_compare,'ALL',all_diff_gene)
        if self.run_kegg:
            kegg_cmd = my_kegg_enrich.run_KEGG_enrich()
            kegg_cmd.extend(my_kegg_enrich.run_kegg_pathview(self.diff_dir, self.diff_program))
            self.cmd_to_scripts(kegg_cmd, self.enrich_thread, my_kegg_enrich.out_dir, self.platform, 'kegg')
        if self.run_go:
            go_cmd = my_go_enrich.run_GO_enrich()
            go_cmd.extend(my_go_enrich.run_GO_DAG(self.topgo))
            self.cmd_to_scripts(go_cmd, self.enrich_thread, my_go_enrich.out_dir, self.platform, 'go')

    def get_mRNA_results(self, root_dir, proj_name):
        ## qc part to be done
        cmd_list = ['#! /bin/sh']
        scripts_dir =os.path.join(root_dir, 'scripts')
        python_tools.circ_mkdir_unix(scripts_dir)
        cp_results_script = os.path.join(scripts_dir, 'cp_results.sh')
        self.result_dir = os.path.join(root_dir, '{0}/mRNA_analysis_results'.format(proj_name))
        cmd_list.append('mkdir -p %s' % self.result_dir)
        diff_gene_list_dir = os.path.join(self.diff_dir,'Diff_list')
        
        qc_result_dir = os.path.join(self.result_dir, 'qc')
        quant_result_dir = os.path.join(self.result_dir, 'quantification')
        diff_result_dir = os.path.join(self.result_dir, 'differential_analysis')
        diff_table_dir = os.path.join(diff_result_dir, 'diff_table')
        diff_vol_dir = os.path.join(diff_result_dir, 'volcano_plot')
        diff_heatmap_dir = os.path.join(diff_result_dir, 'heatmap')
        diff_cluster_dir = os.path.join(diff_result_dir, 'cluster')
        go_result_dir = os.path.join(self.result_dir, 'enrichment_analysis/GO')
        kegg_result_dir = os.path.join(self.result_dir, 'enrichment_analysis/KEGG')
        report_plot_dir = os.path.join(self.result_dir, 'report_plot_tmp')

        ## qc results
        cmd_list.append('echo get qc results start\ndate')
        cmd_list.append('mkdir -p %s' % report_plot_dir)
        cmd_list.append('mkdir -p %s' % qc_result_dir)
        cmd_list.append('cut -f2 %s > %s/sample.list' % (self.group_sample, self.cleandata_dir))
        cmd_list.append('python %s %s/sample.list %s %s' % (FASTQC_SUMMERY, self.cleandata_dir, self.qc_dir, qc_result_dir))
        cmd_list.append('Rscript %s %s/reads_quality %s/reads_quality' % (READS_QUALITY_PLOT, qc_result_dir, qc_result_dir))
        cmd_list.append('echo qc results done\ndate')

        ## quant results
        cmd_list.append('echo cp quantification results start\ndate')
        cmd_list.append('mkdir -p %s' % quant_result_dir)
        cmd_list.append('python %s %s/genes.counts.matrix %s %s/Gene.readcount.xls' % (QUANT_ANNO, self.quant_dir, self.anno_files, quant_result_dir))
        cmd_list.append('python %s --table %s/Gene.readcount.xls' % (DIFF_TABLE_TREAT, quant_result_dir))
        cmd_list.append('cp %s/genes.TMM.EXPR.matrix %s' % (self.quant_dir, quant_result_dir))
        cmd_list.append('python %s --table %s/genes.TMM.EXPR.matrix' % (DIFF_TABLE_TREAT, quant_result_dir))
        cmd_list.append('Rscript %s %s/genes.TMM.EXPR.matrix %s' % (EXP_PLOT, quant_result_dir, quant_result_dir))
        cmd_list.append('Rscript %s %s/genes.TMM.EXPR.matrix %s %s' % (PCA_PLOT, quant_result_dir, self.group_sample, quant_result_dir))
        cmd_list.append('python %s %s/genes.TMM.EXPR.matrix %s %s/Gene.tpm.xls' % (QUANT_ANNO, quant_result_dir, self.anno_files, quant_result_dir))
        cmd_list.append('rm %s/genes.TMM.EXPR.matrix' % quant_result_dir)
        cmd_list.append('cp {self.quant_dir}/sample_correlation/genes.counts.matrix.log2.sample_cor.dat {quant_result_dir}/sample_correlation.data.txt'.format(**locals()))
        cmd_list.append('cp {self.quant_dir}/sample_correlation/genes.counts.matrix.log2.sample_cor_matrix.pdf {quant_result_dir}/sample_correlation.plot.pdf'.format(**locals()))
        cmd_list.append('convert -density 300 {quant_result_dir}/sample_correlation.plot.pdf -quality 90 -flatten {quant_result_dir}/sample_correlation.plot.png'.format(**locals()))
        cmd_list.append('python %s --table %s/sample_correlation.data.txt --add_info Sample' % (DIFF_TABLE_TREAT, quant_result_dir))
#        cmd_list.append('Rscript %s %s/genes.TMM.EXPR.matrix %s' % (EXP_PLOT, self.quant_dir, quant_result_dir))
        cmd_list.append('echo cp quantification results finished\ndate')

        ## cp diff analysis results
        cmd_list.append('echo cp diff analysis results start\ndate')
        for each_compare in self.compare_list:
            each_table_compare_dir = os.path.join(diff_table_dir, each_compare)
            each_vol_compare_dir = os.path.join(diff_vol_dir, each_compare)
            each_de_file = os.path.join(self.diff_dir, 'genes.counts.matrix.%s.%s.DE_results' % (each_compare, self.diff_program))
            each_de_result_file = os.path.join(each_table_compare_dir, '%s.%s.DE_results.txt' % (each_compare, self.diff_program))            
            cmd_list.append('mkdir -p %s' % each_table_compare_dir)            
            cmd_list.append('python %s %s %s %s' % (QUANT_ANNO, each_de_file, self.anno_files, each_de_result_file))
            cmd_list.append('python %s --table %s' % (DIFF_TABLE_TREAT, each_de_result_file))
            cmd_list.append('mkdir -p %s' % each_vol_compare_dir)
            cmd_list.append('Rscript %s %s %s %s %s %s' % (VOLCANO_PLOT, each_de_result_file, each_compare, each_vol_compare_dir, self.qvalue, self.logfc))
            cond_list = each_compare.split('_vs_')            
            for each_cond in cond_list:
                each_name = '%s-UP' % each_cond
                each_cond_up_file = glob.glob(r'%s/*%s.*.DE_results.*.%s-UP.subset' % (self.diff_dir,each_compare,each_cond))[0]
                each_cond_up_result_file = os.path.join(each_table_compare_dir, '%s.%s.subset.txt' % (each_compare, each_name))
                cmd_list.append('python %s %s %s %s' % (QUANT_ANNO, each_cond_up_file, self.anno_files, each_cond_up_result_file))
        heatmap_pdf = os.path.join(self.diff_dir, 'diffExpr.P%s_C%s.matrix.log2.centered.genes_vs_samples_heatmap.pdf' % (self.qvalue, self.logfc))
        cmd_list.append('cp %s/numDE_feature_counts.P%s_C%s.matrix %s/DE_summary.txt' % (self.diff_dir, self.qvalue, self.logfc, diff_table_dir))
        cmd_list.append('python %s --table %s/DE_summary.txt --add_info Group' % (DIFF_TABLE_TREAT, diff_table_dir))        
        cmd_list.append('mkdir -p %s' % diff_heatmap_dir)
        cmd_list.append('cp %s %s/Heatmap.pdf' % (heatmap_pdf, diff_heatmap_dir))
        cmd_list.append('convert -density 300 {diff_heatmap_dir}/Heatmap.pdf -quality 90 -flatten {diff_heatmap_dir}/Heatmap.png'.format(**locals()))
        cmd_list.append('python %s %s/diffExpr.P%s_C%s.matrix %s %s/diff_gene.tpm.txt' % (QUANT_ANNO, self.diff_dir, self.qvalue, self.logfc, self.anno_files, diff_heatmap_dir))
        cmd_list.append('python %s --table %s/diff_gene.tpm.txt' % (DIFF_TABLE_TREAT, diff_heatmap_dir))                
        cmd_list.append('mkdir -p %s' % diff_cluster_dir)
        for per in [10, 15, 20]:
            cmd_list.append('cp -r %s/diffExpr.P%s_C%s.matrix.RData.clusters_fixed_P_%s %s/clusters_fixed_P_%s' % (self.diff_dir, self.qvalue, self.logfc, per, diff_cluster_dir, per))
            cmd_list.append('rm %s/clusters_fixed_P_%s/*R' % (diff_cluster_dir, per))
        cmd_list.append('Rscript %s %s/clusters_fixed_P_20 %s' % (CLUSTER_LINE, diff_cluster_dir, report_plot_dir))
        cmd_list.append('echo cp diff analysis results finished\ndate')

        ## cp go enrichment analysis results
        cmd_list.append('echo cp go enrichment analysis results start\ndate')
        for each_compare in self.compare_list:
            ## cp bar plot
            each_go_bar_dir = os.path.join(go_result_dir, 'GO_bar_plot/%s' % each_compare)
            cmd_list.append('mkdir -p %s' % (each_go_bar_dir))
            cmd_list.append('cp %s/GO/%s/*bar* %s' % (self.enrich_dir, each_compare, each_go_bar_dir))
            ## cp dag plot
            each_go_dag_dir = os.path.join(go_result_dir, 'GO_dag_plot/%s' % each_compare)
            cmd_list.append('mkdir -p %s' % (each_go_dag_dir))
            cmd_list.append('cp -r %s/GO/%s/GO_DAG/* %s' % (self.enrich_dir, each_compare, each_go_dag_dir))
            ## cp enrichment table
            each_go_tabel_dir = os.path.join(go_result_dir, 'GO_enrich_table/%s' % each_compare)
            cmd_list.append('mkdir -p %s' % (each_go_tabel_dir))
            cond_list = ['ALL']
            cond_list.extend(['%s-UP' % each for each in each_compare.split('_vs_')])
            for each_cond in cond_list:
                each_go_file = os.path.join(self.enrich_dir, '%s/GO/%s/%s.%s.GO.enrich.xls' % (self.enrich_dir, each_compare, each_compare, each_cond))
                each_go_out = os.path.join(self.enrich_dir, '%s/%s.%s.GO.enrich.xls' % (each_go_tabel_dir, each_compare, each_cond))
                diff_gene_list_dir = os.path.join(self.diff_dir,'Diff_list')
                each_diff_list = os.path.join(diff_gene_list_dir, '%s.%s.list' % (each_compare, each_cond))
                cmd_list.append('python %s %s %s %s %s' % (GO_RESULT_ANNO, each_go_file, each_diff_list, self.go, each_go_out))
        cmd_list.append('echo cp go enrichment analysis results finished\ndate')

        ## cp kegg enrichment analysis results
        cmd_list.append('echo cp kegg enrichment analysis results start\ndate')
        for each_compare in self.compare_list:
            ## cp kegg bar plot
            each_kegg_bar_dir = os.path.join(kegg_result_dir, 'KEGG_bar_plot/%s' % each_compare)
            cmd_list.append('mkdir -p %s' % (each_kegg_bar_dir))
            cmd_list.append('cp %s/KEGG/%s/*bar* %s' % (self.enrich_dir, each_compare, each_kegg_bar_dir))
            ## cp kegg table
            each_kegg_tabel_dir = os.path.join(kegg_result_dir, 'KEGG_enrich_table/%s' % each_compare)
            cmd_list.append('mkdir -p %s' % (each_kegg_tabel_dir))
            cmd_list.append('cp %s/KEGG/%s/*enrich.xls %s' % (self.enrich_dir, each_compare, each_kegg_tabel_dir))
            ## cp kegg pathway
            each_kegg_pathway_dir = os.path.join(kegg_result_dir, 'KEGG_pathway/%s' % each_compare)
            cmd_list.append('mkdir -p %s' % (each_kegg_pathway_dir))
            cmd_list.append('cp -r %s/KEGG/%s/*pathway %s' % (self.enrich_dir, each_compare, each_kegg_pathway_dir))
        cmd_list.append('echo cp kegg enrichment analysis results finished\ndate')
        ## finished
        python_tools.write_obj_to_file(cmd_list, cp_results_script)

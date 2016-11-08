'''
data:2016-3-26

RNAseq tools and functions

last modification:2016-3-26
'''

import sys
import os
import python_tools
import gzip
from HTSeq import GFF_Reader

## database
KO_PEP_DIR = '/home/lxgui/test_kobas/seq_pep/'

## softwares and scripts absolute path
R_BIN = '/usr/local/bin/'
R_BIN_3_1 = '/home/public/R_packages/R-3.1.0/bin/'
BLAST_BIN = '/usr/bin/'
KALLISTO = '/home/public/software/kallisto-0.42.4/kallisto'


GOSEQ = '/home/lxgui/scripts/goseq.R'
ENRICH_BAR = '/home/lxgui/scripts/enrich_bar.2016-3-31.R'

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

def get_gene_info(transcript_info_dict) :
    gene_info_dict = {}
    for each_tr in transcript_info_dict :
        gene_id = transcript_info_dict[each_tr]
        if gene_id in gene_info_dict :
            each_tr_len = transcript_info_dict[each_tr]['length']
            each_tr_exon = transcript_info_dict[each_tr]['exon_num']
            each_tr_tss = transcript_info_dict[each_tr]['start']
            each_tr_tts = transcript_info_dict[each_tr]['end']
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
    return gene_info_dict

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
        self.output = None
        self.auot_run = True
        self.plot_data = None
        self.compare_name = None

    def get_target_length_table(self):
        if os.path.isfile(self.target_length) :
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
                    each_gene_len = python_tools.Median(self.transcript_dict[each_gene]['transcript_len'])
                    target_length_info.write('%s\t%s\n' % (each_gene,each_gene_len))
            target_length_info.close()

    def run_goseq_enrich(self) :
        cmd = '%s/Rscript %s %s %s %s %s' % (R_BIN,GOSEQ,self.target_list,self.target_length,self.go,self.output)
        if self.auot_run :
            python_tools.circ_call_process(cmd)
        return cmd

    def run_plot(self) :
        out_dir = os.path.split(self.output)[0]
        cmd = '%s/Rscript %s %s %s %s %s' % (R_BIN_3_1,ENRICH_BAR,self.compare_name,self.output,'GO',out_dir)
        if self.auot_run :
            python_tools.circ_call_process(cmd)
        return cmd

class KEGG_enrich :
    def __init__(self) :
        self.seq = None
        self.all_blast_out = None
        self.compare_blast_out = None
        self.species = None        
        self.ko_seq = None  
        self.auot_run = True
        self.plot_data = None
        self.compare_name = None
        self.blast_dict = {}
        self.target_type = 'gene'
        self.gtf = None
        self.transcript_dict = {} 
        self.diff_list = None
        self.output = None
        self.out_dir = None
        self.database = 'Ensembl'

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
        if self.auot_run :
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

        # diff_blast = os.path.join(out_dir,'%s.blasttab.json' % species)
        # with open(diff_blast,'w') as diff_blast_file :
        #     json.dump(self.blast_dict,diff_blast_file)
        # return diff_blast     

    def extract_diff_blast(self):
        # with open(gene_blast_json,'r') as gene_blast_info :
        #     gene_blast_dict = json.load(gene_blast_info)
        self.out_dir = os.path.split(self.output)[0]
        self.compare_blast_out = os.path.join(self.out_dir,'%s.blasttab' % self.compare_name)
        with open(self.compare_blast_out,'w') as diff_blast_out :
            for each_gene in self.diff_list :
                if each_gene in self.blast_dict :
                    blast_info = '\t'.join(self.blast_dict[each_gene])
                    diff_blast_out.write('%s\t%s\n' % (each_gene,blast_info))

    def KEGG_enrich(self):        
        if self.database == 'NCBI' :
            cmd = 'run_kobas.py -i {self.diff_list}  -t id:ncbigene -s {self.species} -d K -o {self.output}'.format(**locals())
        else :
            cmd = 'run_kobas.py -i {self.compare_blast_out}  -t blastout:tab -s {self.species} -d K -o {self.output}'.format(**locals())
        if self.auot_run :
            python_tools.circ_call_process(cmd)
        return cmd    

    def treat_KEGG_table(self):
        kegg_out_dir,kegg_out_name = os.path.split(self.output)
        kegg_tmp_file = os.path.join(kegg_out_dir,'tmp.%s' % kegg_out_name)
        os.system('mv %s %s' % (self.output,kegg_tmp_file))
        kegg_out_info = open(self.output,'w')
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

    def run_plot(self) :
        cmd = '%s/Rscript %s %s %s %s %s' % (R_BIN_3_1,ENRICH_BAR,self.compare_name,self.output,'KEGG',self.out_dir)
        if self.auot_run :
            python_tools.circ_call_process(cmd)
        return cmd

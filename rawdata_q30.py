#coding=utf-8
'''

calculate q30 of all rawdata for a given project/key word

'''

import argparse
from os.path import join as join_path
import xlsxwriter
from os import getcwd as getcwd

cwd = getcwd()
all_seq_info_file = '/home/lxgui//Documents/paiji/check_up_dowm/出库样品信息.txt'
seq_info_file = '/home/lxgui//Documents/paiji/proj_q30_summary/下机整合.txt'

parser = argparse.ArgumentParser()
parser.add_argument('--seq_info', help = 'sequencing results table', default = seq_info_file)
parser.add_argument('--proj_name', help = 'project names, sep with ","', required = True)
parser.add_argument('--all_seq_info', help = 'all samples sequencing info', default = all_seq_info_file)
parser.add_argument('--seq_num', help = 'seq numbers, sep with ","', default = "")
parser.add_argument('--out_dir', help = 'project rawdata directory', required = True)
parser.add_argument('--ALL', help = 'sample info in all seq numbers ', action='store_true', default = False)
args = parser.parse_args()

seq_list = args.seq_num.split(',')
proj_name_list =  args.proj_name.split(',')

def output_to_xlsx(worksheet,info_list,row_num,format  = None):    
    for n,each in enumerate(info_list):
        each = str(each)
        if format:
            worksheet.write_string(row_num,n,each,format)
        else:
            worksheet.write_string(row_num,n,each)

sample_seq_dict = {}
with open(args.seq_info) as seq_info_file :
    for n,eachline in enumerate(seq_info_file) :
        if n != 0 :
            eachline_info = eachline.strip().split('\t')
            wgc_id = eachline_info[1].strip()
            q30 = float(eachline_info[-2])
            seq_data = float(eachline_info[-4])
            seq_num = eachline_info[-1]
            reads_num = float(eachline_info[-3])
            if not args.ALL and seq_num not in seq_list:
                continue
            if wgc_id not in sample_seq_dict :
                sample_seq_dict[wgc_id] = [[seq_data],[q30],[reads_num]]
            else :
                sample_seq_dict[wgc_id][0].append(seq_data)
                sample_seq_dict[wgc_id][1].append(q30)
                sample_seq_dict[wgc_id][2].append(reads_num)

output_file = join_path(args.out_dir,'Data_summary.xlsx')
output_file2 = join_path(args.out_dir,'Data_summary.txt')
output_file2_info =  open(output_file2,'w')

workbook = xlsxwriter.Workbook(output_file)
bold = workbook.add_format({'bold': 1})
worksheet = workbook.add_worksheet()
with open(args.all_seq_info) as proj_file :
    row_num = 0
    for n,eachline in enumerate(proj_file) :
        eachline = eachline.strip()        
        if n == 0 :
            output_file2_info.write('Sample_ID\tIndex\tReads_number(M)\tData_size(G)\tLength\tQ30(%)\n')
            output_info_list = ['Sample_ID','Index','Reads_number(M)','Data_size(G)','Length','Q30(%)']
            output_to_xlsx(worksheet,output_info_list,row_num,bold)
            row_num += 1
        else :
            eachline_info = eachline.split('\t')
            each_proj_name = eachline_info[0]
            if each_proj_name in proj_name_list:
                wgc_id = eachline_info[1]
                if wgc_id not in sample_seq_dict:
                    continue
                sample_id = eachline_info[4]
                sample_index = eachline_info[7]
                old_data_size = float(eachline_info[15])
                data_size = sum(sample_seq_dict[wgc_id][0])
                output_data_size = round(data_size,2)
                reads_num_sum = round(sum(sample_seq_dict[wgc_id][2]),2)
                if output_data_size != round(old_data_size,2) :
                    print '%s different data_size : %s(now),%s(before)' % (wgc_id,data_size,old_data_size)
                q30_sum = 0
                for n,each_q30 in enumerate(sample_seq_dict[wgc_id][1]) :
                    q30_sum += each_q30 * sample_seq_dict[wgc_id][0][n]
                q30 = round(q30_sum/data_size,2)                
                # output_info.write('%s\t%s\t%s\t%s\tPE150\t%s\n' % (sample_id,sample_index,reads_num_sum,data_size,q30))
                output_info_list = [sample_id,sample_index,reads_num_sum,output_data_size,'PE150',q30]
                output_to_xlsx(worksheet,output_info_list,row_num)
                row_num += 1
                output_file2_info.write('{sample_id}\t{sample_index}\t{reads_num_sum}\t{output_data_size}\tPE150\t{q30}\n'.format(**locals()))
workbook.close()
output_file2_info.close()

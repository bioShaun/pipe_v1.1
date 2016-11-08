#coding=utf-8

'''
generate sequencing report (pdf)

'''
import argparse
import os
import sys
reload(sys)
sys.setdefaultencoding('utf-8')

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--project_name", required = True)
parser.add_argument("-ss", "--sample_size", default = '')
parser.add_argument("-ds", "--data_size_sample")
parser.add_argument("-da", "--data_size_all")
parser.add_argument("-an", "--add_sample_num", default = '')
parser.add_argument("-t", "--data_table", required = True)
parser.add_argument('-o','--out_dir', required = True)
args = parser.parse_args()

if os.path.exists(args.out_dir):
    os.system('mkdir -p %s' % args.out_dir)

tb_info = open(args.data_table).readlines()
if not args.sample_size:
    all_sample_size = len(tb_info)
else:
    all_sample_size = int(args.sample_size)
seq_sample_size = len(tb_info)

template = '/home/lxgui/Documents/latex_test/report_pdf/sequencing_report_template.tex'
output_file = os.path.join(args.out_dir,'%s.tex' % args.project_name)

output_file_info = open(output_file,'w')
with open(template) as template_info:
    for eachline in template_info:
        space = eachline.split('%')[0]
        eachline_test = eachline.lstrip()
        if eachline_test.startswith('%project_id'):
            output = '%s\\textbf{项目编号:} %s \\\\[1cm]\n' % (space,args.project_name)
        elif eachline_test.startswith('%project_summery'):            
            output = '''
%s\\item \\textbf{项目任务步骤:}X Ten模式，150PE，
完成共计%s个文库建库及上机，
''' % (space, all_sample_size)
            if args.data_size_sample:
                output += '%s每个文库产出数据量%sG。\n' % (space, args.data_size_sample)
            elif args.data_size_all:
                output += '%s共产出数据量%sG。\n' % (space, args.data_size_all)
            else:
                sys.exit('--data_size_sample or --data_size_all is needed!')
        elif eachline_test.startswith('%project_info'):
            output = '%s\\item \\textbf{完成情况:} 已经完成%s个文库。\n' % (space, seq_sample_size)
        elif eachline_test.startswith('%project_sup') and args.add_sample_num:
            output = '%s\\item \\textbf{说明:} 7个文库正在补测。\n' % (space, seq_sample_size)
        elif eachline_test.startswith('%table_example'):
            table_info_list = []
            for each_row in tb_info:
                each_row = each_row.strip().split('\t')
                each_row_info = ' & '.join(each_row)
                each_row_info = '%s%s\\\\' % (space, each_row_info)
                table_info_list.append(each_row_info)
                table_info_list.append('%s\\hline' % (space))
            output = '%s\n' % '\n'.join(table_info_list) 
        else:
            output = eachline
        output_file_info.write(output)
output_file_info.close()



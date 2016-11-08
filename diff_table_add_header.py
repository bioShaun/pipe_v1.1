import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--table', help = 'differentail analysis result', required = True)
parser.add_argument('--add_info', help = 'info need to add in table header', default = 'id')
args = parser.parse_args()

tmp_file = '%s.bak' % args.table

os.system('mv %s %s' % (args.table,tmp_file))

new_table = open(args.table,'w')
with open(tmp_file) as tmp_file_info:
    for n,eachline in enumerate(tmp_file_info):
        if n == 0 and (not eachline.startswith(args.add_info)):
            eachline = eachline.lstrip()
            new_table.write('%s\t%s' % (args.add_info,eachline))
        else:
            new_table.write(eachline)
new_table.close()
os.system('rm %s' % (tmp_file))


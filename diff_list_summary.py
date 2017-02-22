'''
diff gene summary
column1: gene_ID
column2: differential expressed in which compare

'''

import sys
import os
import glob
import python_tools
import re

if not len(sys.argv) == 2:
    print '    python ' + sys.argv[0] + ' diff.list.dir '
    sys.exit(0)

diff_list_dir = sys.argv[1]

diff_list_files = glob.glob('{}/*UP.list'.format(diff_list_dir))

def get_compare_staus(diff_list_name):
    compare_name1, tmp = diff_list_name.split('_vs_')
    if compare_name1 in tmp:
        compare_name2 = re.match('(\S+).{}-UP.list'.format(compare_name1),tmp).groups()[0]
        return '{compare_name1}_vs_{compare_name2}'.format(**locals()),'{compare_name1}-UP'.format(**locals())
    else:
        compare_name2_twice = re.match('(\S+)-UP.list',tmp).groups()[0]
        compare_name2 = compare_name2_twice[0:(len(compare_name2_twice)-1)/2]
        return '{compare_name1}_vs_{compare_name2}'.format(**locals()),'{compare_name2}-UP'.format(**locals())

diff_gene_dict = {}
for each_file in diff_list_files:
    each_file_name = os.path.basename(each_file)
    compare, detail = get_compare_staus(each_file_name)
    each_file_gene_list = [each.strip() for each in open(each_file)]
    for each_gene in each_file_gene_list:
        if each_gene in diff_gene_dict:
            diff_gene_dict[each_gene].append('{compare}({detail})'.format(**locals()))
        else:
            diff_gene_dict[each_gene] = ['{compare}({detail})'.format(**locals())]

diff_gene_summary_file = os.path.join(diff_list_dir, 'diffgene.summary.txt')
with open(diff_gene_summary_file, 'w') as diff_gene_summary_file_inf:
    diff_gene_summary_file_inf.write('Gene_ID\tDifferential_expressed_groups\n')
    for each_gene in diff_gene_dict:
        each_diff_stat = ','.join(diff_gene_dict[each_gene])
        diff_gene_summary_file_inf.write('{0}\t{1}\n'.format(each_gene, each_diff_stat))


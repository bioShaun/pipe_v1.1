import sys
import os
import python_tools
import glob
import subprocess

if not len(sys.argv) == 5:
    print '    python ' + sys.argv[0] + ' diff.analysis.dir compare.list gene.anno diff.result.dir'
    sys.exit(0)

diff_dir = sys.argv[1]
compare_file = sys.argv[2]
gene_anno = sys.argv[3]
out_dir = sys.argv[4]

compare_list = [each.strip() for each in open(compare_file)]

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
    ret_code = p.wait()
    output = p.communicate()[0]
    return output

## cp diff table
for each in compare_list:
    each_de_results = os.path.join(diff_dir, 'genes.counts.matrix.{0}.edgeR.DE_results'.format(each))
    each_de_results_dir = os.path.join(out_dir, each)
    if not os.path.exists(each_de_results_dir):
        python_tools.circ_mkdir_unix(each_de_results_dir)
    each_de_results_out_tmp = os.path.join(each_de_results_dir, 'tmp.{0}.edgeR.DE_results.txt'.format(each))
    each_de_results_out = os.path.join(each_de_results_dir, '{0}.edgeR.DE_results.txt'.format(each))
    ## add table header
    cmd_list  = []
    cmd_list.append(['cp', each_de_results, each_de_results_out_tmp])
    cmd_list.append(['python', '/home/lxgui/scripts/diff_table_add_header.py', '--table', each_de_results_out_tmp, '--add_info', 'Gene_ID'])
    ## vocalno plot
    cmd_list.append(['Rscript', '/home/lxgui/scripts/Volcano_Plot_20160406.R', each_de_results_out_tmp, each, each_de_results_dir, '0.001', '2'])
    cmd_list.append(['python', '/home/lxgui/scripts/add_gene_anno_v2.py', each_de_results_out_tmp, gene_anno, each_de_results_out])
    cmd_list.append(['rm', each_de_results_out_tmp])
    each_sub_reg_list = each.split('_vs_')
    for each_sub in each_sub_reg_list:
        name = '%s-UP' % each_sub
        each_sub_diff_result = glob.glob(r'{0}/genes.counts.matrix.{1}.edgeR.DE_results.*.{2}.subset'.format(diff_dir, each, name))[0]
        each_sub_diff_out = os.path.join(each_de_results_dir, '{0}.{1}.subset.txt'.format(each, name))
        cmd_list.append(['python', '/home/lxgui/scripts/add_gene_anno_v2.py', each_sub_diff_result, gene_anno, each_sub_diff_out])
    for each_cmd in cmd_list:
        log_tmp = run_cmd(each_cmd)

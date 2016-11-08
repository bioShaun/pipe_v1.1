import os
import sys
import python_tools

if not len(sys.argv) == 5:
    print '    python ' +sys.argv[0] + ' compare_list_file enrich_dir go_file out_dir'
    sys.exit(0)

compare_list_file = sys.argv[1]
enrich_dir = sys.argv[2]
go_file = sys.argv[3]
out_dir = sys.argv[4]

compare_list = [each.strip() for each in open(compare_list_file)]

for each in compare_list:
    each_enrich = os.path.join(enrich_dir, each)
    each_enrich_out = os.path.join(out_dir, each)
    if not os.path.exists(each_enrich_out):
        python_tools.circ_mkdir_unix(each_enrich_out)
    reg_list = ['ALL']
    each_sub_reg_list = each.split('_vs_')
    reg_list.extend(each_sub_reg_list)
    for each_sub in reg_list:
        name = each_sub
        if each_sub != 'ALL':
            name = '%s-UP' % each_sub
        each_diff_file1 = os.path.join(each_enrich, '%s-target.list' % name) 
        each_diff_file2 = os.path.join(each_enrich, '%s.list' % name)
        if os.path.exists(each_diff_file1):
            each_diff_file = each_diff_file1
        else:
            each_diff_file = each_diff_file2
        each_go_file1 = os.path.join(each_enrich, '%s.%s.GO.enrich.xls' % (each, name))
        each_go_file2 = os.path.join(each_enrich, '%s.%s-target.GO.enrich.xls' % (each, name))
        if os.path.exists(each_go_file2):
            each_go_file = each_go_file2
            each_go_annnotate = os.path.join(each_enrich_out, '%s.%s-target.GO.enrich.xls' % (each, name))
        else:
            each_go_file = each_go_file1
            each_go_annnotate = os.path.join(each_enrich_out, '%s.%s.GO.enrich.xls' % (each, name))
        print ('python ~/scripts/GO_add_diff_gene.py %s %s %s %s' % (each_go_file, each_diff_file, go_file, each_go_annnotate))

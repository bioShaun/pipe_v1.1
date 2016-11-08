import sys
import os

if not len(sys.argv) == 4:
    print 'python ' + sys.argv[0] + ' samples expression out_dir'
    sys.exit(0)

sample_list_file = sys.argv[1]
exp_file = sys.argv[2]
out_dir = sys.argv[3]

sample_list = [each.strip() for each in open(sample_list_file)]

out_put_gene_file = os.path.join(out_dir, 'expressed.id.list')
out_put_exp_file = os.path.join(out_dir, 'expressed.table.txt')

sample_index = []
header = ['ID']
out_put_gene_file_info = open(out_put_gene_file, 'w')
out_put_exp_file_info = open(out_put_exp_file, 'w')
with open(exp_file) as exp_file_info:
    for n, eachline in enumerate(exp_file_info):
        eachline_info = eachline.strip().split('\t')
        if n == 0:
            for m, each in enumerate(eachline_info):
                if each in sample_list:
                    sample_index.append(m)
                    header.append(each)
            out_put_exp_file_info.write('%s\n' % '\t'.join(header))
        else:
            each_id = eachline_info[0]
            exp_list = []
            for m, each in enumerate(eachline_info):
                if m in sample_index:
                    exp_list.append(float(each))
            if max(exp_list) > 0.1:
                exp_list = [str(each) for each in exp_list]
                exp_out = '\t'.join(exp_list)
                out_put_exp_file_info.write('%s\t%s\n' % (each_id, exp_out))
                out_put_gene_file_info.write('%s\n' % each_id)
out_put_gene_file_info.close()
out_put_exp_file_info.close()


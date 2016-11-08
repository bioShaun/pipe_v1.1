import sys
import os

if not len(sys.argv) == 3:
    print '    python ' + sys.argv[0] + ' exp.table out.dir'
    sys.exit(0)

def write_list_to_file(obj, fn, header = ''):
    fh = open(fn, 'w')
    if header != '':
        fh.write('%s\n' % header)
    for each_list in obj:
        each_list = [str(each) for each in each_list]
        eachline = '\t'.join(each_list)
        fh.write('%s\n' % eachline)
    fh.close()  

exp_table = sys.argv[1]
out_dir = sys.argv[2]

inter_list = [0.5, 1, 2, 4, 8]
header = ''
out_put_list = [[], [], [], [], []]
with open(exp_table) as exp_table_inf:
    for n, eachline in enumerate(exp_table_inf):
        if n == 0:
            header = eachline.strip()
        else:
            eachline_inf = eachline.strip().split('\t')
            eachline_exp = [ float(each) for each in eachline_inf[1:] ]
            max_exp = max(eachline_exp)
            for m, each in enumerate(inter_list):
                if max_exp <= each:
                    if m == 0:
                        break
                    else:
                        out_put_list[m-1].append(eachline_inf)
                        break
            else:
                out_put_list[-1].append(eachline_inf)

for n, each in enumerate(inter_list):
    if n+1 < len(inter_list):
        out_name = 'tpm_between_%s_%s.txt' % (inter_list[n], inter_list[n+1])
    else:
        out_name = 'tpm_gt_%s.txt' % inter_list[n]
    out_file = os.path.join(out_dir, out_name)
    write_list_to_file(out_put_list[n], out_file, header)

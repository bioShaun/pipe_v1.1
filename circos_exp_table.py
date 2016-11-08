import sys
import os
import python_tools

if not len(sys.argv) == 5:
    print '    python ' + sys.argv[0] + ' chr.list bin_width exp.table circos.exp.table'
    sys.exit(0)

def list_add(list1,list2):
    list_sum = []
    for n, each in enumerate(list1):
        list_sum.append(list1[n] + list2[n])
    return list_sum

def write_tuple_to_file(fn, chrom, bin_width, chr_tuple, chr_dict, append = True):
    fh = open(fn, 'a' if append is True else 'w')
    chrom = str(chrom)
    output_chr = 'chr%s' % chrom
    chr_end = int(chr_dict[chrom])
    for each in chr_tuple:
        start = each[0]*bin_width
        end = (each[0]+1)*bin_width
        if end > chr_end:
            end = chr_end
        exp_info = [str(each) for each in each[1]]
        exp_out = '\t'.join(exp_info)
        fh.write('%s\t%s\t%s\t%s\n' % (output_chr, start, end, exp_out))
    fh.close()

chr_list_file = sys.argv[1]
bin_width = int(sys.argv[2])
exp_table = sys.argv[3]
output = sys.argv[4]

chr_dict = python_tools.table_to_dict(chr_list_file, 1, 2, False)

chr_bin_exp_dict = {}

chr_tmp = ""
with open(exp_table) as exp_table_info:
    for n, eachline in enumerate(exp_table_info):
        if n != 0:
            eachline_info = eachline.strip().split('\t')
            chrom = eachline_info[0]
            if chr_tmp == "":
                chr_tmp = chrom
            if str(chrom) not in chr_dict:
                continue
            start = int(eachline_info[1])
            end = int(eachline_info[2])
            exp_data = [float(each) for each in eachline_info[5:]]
            chr_index1 = start/bin_width
            chr_index2 = end/bin_width
            if chrom != chr_tmp :
                sorted_chr_bin_exp_dict = sorted(chr_bin_exp_dict.iteritems(), key=lambda x : x[0])
                write_tuple_to_file(output, chr_tmp, bin_width, sorted_chr_bin_exp_dict, chr_dict)
                chr_bin_exp_dict = {}
                chr_tmp = chrom           
            for m in range(chr_index1,chr_index2+1):
                if m not in chr_bin_exp_dict:
                    chr_bin_exp_dict[m] = exp_data
                else:
                    chr_bin_exp_dict[m] = list_add(chr_bin_exp_dict[m], exp_data)

sorted_chr_bin_exp_dict = sorted(chr_bin_exp_dict.iteritems(), key=lambda x : x[0])
write_tuple_to_file(output, chr_tmp, bin_width, sorted_chr_bin_exp_dict, chr_dict)
chr_bin_exp_dict = {}



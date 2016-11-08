import os
import sys

if not len(sys.argv) == 3:
    print 'python ' + sys.argv[0] + ' file_dir summary'
    sys.exit(0)

file_dir = sys.argv[1]
summary = sys.argv[2]

file_list = os.listdir(file_dir)
summary_info = open(summary, 'w')
summary_info.write('File_name\tIndex\tAllele1\tAllele1_avg\tAllele2\tAllele2_avg\tMean_Sq\tF-value\tP-value\tAllele1_count\tAllele2_count\tSignificant\n')

sort_file_list = sorted(file_list)

for each_file in sort_file_list:
    each_file_path = os.path.join(file_dir, each_file)
    each_file_name = '_'.join(each_file.split('_')[0:2])
    each_file_info1_list = []
    each_file_info2_list = []
    allele1 = ''
    allele2 = ''
    allele1_index = 0
    allele2_index = 0
    with open(each_file_path) as each_file_info:
        sig = 'No'
        for n, eachline in enumerate(each_file_info):            
            eachline_info = eachline.strip().split()
            if eachline.startswith('Signif'):
                sig = 'Yes'
            if sig == 'Yes':
                m = n-2
            else:
                m = n
            if m == 0:
                snp_id = eachline.strip().split('`')[1]
                each_file_info1_list.append(snp_id)
            elif m == 2:
                mean = eachline_info[3]
                F_value = eachline_info[4]
                Pvalue = eachline_info[5]
                each_file_info2_list.extend([mean, F_value, Pvalue])
            elif m == 7:
                allele1 = eachline_info[1]
                allele1_avg = eachline_info[2]
                each_file_info1_list.extend([allele1, allele1_avg])
            elif m == 8:
                allele2 = eachline_info[1]
                allele2_avg = eachline_info[2]
                each_file_info1_list.extend([allele2, allele2_avg])
            elif m == 12:
                allele1_index = eachline_info.index(allele1)
                allele2_index = eachline_info.index(allele2)
            elif m == 13:
                allele1_count = eachline_info[allele1_index]
                allele2_count = eachline_info[allele2_index]
                each_file_info1_list.extend([allele1_count, allele2_count])
    eachfile_out = '%s\t%s\t%s\t%s\n' % (each_file_name,'\t'.join(each_file_info1_list), '\t'.join(each_file_info2_list), sig)
    summary_info.write(eachfile_out)
summary_info.close()


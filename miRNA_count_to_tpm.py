import sys

if not len(sys.argv) == 3:
    print 'python ' + sys.argv[0] + ' miRNA_count_table miRNA_tpm_table'
    sys.exit(0)

miRNA_count_table = sys.argv[1]
miRNA_tpm_table = sys.argv[2]

miRNA_tpm_table_info = open(miRNA_tpm_table, 'w')
miRNA_count_dict = {}
gene_list = []
row_num_list = []

with open(miRNA_count_table) as miRNA_count_table_info:
    for n, eachline in enumerate(miRNA_count_table_info):
        eachline_info = eachline.strip().split('\t')
        if n == 0:
            miRNA_tpm_table_info.write(eachline)
        else:
            for m, each in enumerate(eachline_info):
                if m == 0:
                    gene_list.append(each)
                else:
                    if m not in row_num_list:
                        row_num_list.append(m)
                    if m not in miRNA_count_dict:
                        miRNA_count_dict[m] = [int(each)]
                    else:
                        miRNA_count_dict[m].append(int(each))

for n, each_gene in enumerate(gene_list):
    each_tpm_list = []
    for m, each_row in enumerate(row_num_list):
        row_num = m+1
        each_tpm = str(10**6*miRNA_count_dict[row_num][n]/float(sum(miRNA_count_dict[row_num])))
        each_tpm_list.append(each_tpm)
    miRNA_tpm_table_info.write('%s\t%s\n' % (each_gene, '\t'.join(each_tpm_list)))
miRNA_tpm_table_info.close()


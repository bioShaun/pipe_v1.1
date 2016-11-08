import sys
import os

if not len(sys.argv) == 4:
    print '    python ' + sys.argv[0] + ' miRNA.target miRNA.cor miRNA.tar.cor'
    sys.exit(0)

miRNA_tar = sys.argv[1]
miRNA_cor = sys.argv[2]
miRNA_out = sys.argv[3]

miRNA_pair_cor_dict = {}
with open(miRNA_cor) as miRNA_cor_inf:
    for eachline in miRNA_cor_inf:
        eachline_inf = eachline.strip().split('\t')
        each_miRNA = eachline_inf[0]
        each_mRNA = eachline_inf[1]
        each_cor = eachline_inf[2]
        each_pval = eachline_inf[3]
        miRNA_pair_cor_dict.setdefault(each_miRNA, {})[each_mRNA] = [each_cor, each_pval]

miRNA_out_inf = open(miRNA_out, 'w')
with open(miRNA_tar) as miRNA_tar_inf:
    for n, eachline in enumerate(miRNA_tar_inf):
        eachline = eachline.strip()
        if n == 0:
            miRNA_out_inf.write('%s\tPCC\tPvalue(PCC)\n' % eachline)
        else:
            eachline_inf = eachline.strip().split('\t')
            cor, pval = miRNA_pair_cor_dict[eachline_inf[0]][eachline_inf[1]]
            miRNA_out_inf.write('%s\t%s\t%s\n' % (eachline, cor, pval))
miRNA_out_inf.close()

import sys
import os

if len(sys.argv) != 3:
    print 'python '+ sys.argv[0] + ' mirdeep2_novel_prediction output_dir'
    sys.exit(0)

novel_miRNA_seq_file = sys.argv[1]
output_dir = sys.argv[2]

mature_fa = os.path.join(output_dir, 'novel.mature.fa')
star_fa = os.path.join(output_dir, 'novel.star.fa')
hairpin_fa = os.path.join(output_dir, 'novel.hairpin.fa')
seq_file_list = [mature_fa, star_fa, hairpin_fa]
seq_file_info_list = [open(each,'w') for each in seq_file_list]

with open(novel_miRNA_seq_file) as novel_miRNA_seq_file_info:
    flag = 0
    for n, eachline in enumerate(novel_miRNA_seq_file_info):                
        if 'novel miRNAs predicted by miRDeep2' in eachline:
            flag = 1
        if 'mature miRBase miRNAs detected by miRDeep2' in eachline:
            flag = 0
        eachline_info = eachline.strip().split('\t')
        if flag == 1 and len(eachline_info) == 17:
            if eachline_info[0] == 'provisional id':
                continue
            novel_id = "novel:%s" % eachline_info[0]
            seq_list = eachline_info[13:16]
            for n, each_seq in enumerate(seq_list):
                seq_file_info_list[n].write('>%s\n%s\n' % (novel_id, each_seq))

[each.close() for each in seq_file_info_list]


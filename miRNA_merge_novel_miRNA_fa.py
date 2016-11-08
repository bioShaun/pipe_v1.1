import argparse
import os
import sys
import glob
import re
import python_tools

parser = argparse.ArgumentParser()
parser.add_argument('-in', '--novel_miRNA_file', dest="novel_miRNA_file", required = True)
parser.add_argument('-d', '--mirdeep2_dir', dest = 'mirdeep2_dir', required = True)
parser.add_argument('-o', '--output_dir', dest = 'output_dir', required = True)
args = parser.parse_args()

def merge_overlap_seq(a = 'ggaccaggcuucauuccc', b = 'accaggcuucauucccaa', step = 11):

    '''
    merge overlaped miRNA sequence

    '''

    merged_seq = ''
    if len(a) >= len(b):
        long_seq = a
        short_seq = b
    else:
        long_seq = b
        short_seq = a

    if short_seq in long_seq:
        merged_seq = long_seq
    elif len(short_seq) > step:
        if short_seq[0:step] not in long_seq:
            for m in range(1,len(short_seq)- step + 2):
                if short_seq[m:m + step] in long_seq:
                    header = ''.join(short_seq[0:m])
                    merged_seq = '%s%s' % (header, long_seq)
                    break
            else:
                sys.exit('"%s" and "%s" overlap less than 11bp' % (long_seq, short_seq))
                merged_seq = long_seq
        else:
            for m in range(1, len(short_seq)- step + 2):
                if short_seq[-step - m: -m] in long_seq:
                    tailer = ''.join(short_seq[-m:])
                    merged_seq = '%s%s' % (long_seq, tailer)
                    break
            else:
                sys.exit('"%s" and "%s" overlap less than 11bp' % (long_seq, short_seq))
                merged_seq = long_seq
    else:
        sys.exit('"%s" and "%s" overlap less than 11bp' % (long_seq, short_seq))
        merged_seq = long_seq

    return merged_seq

def merge_seq_list(seq_list):
    merged_seq = ''
    seq_list_cp = seq_list[:]
    if len(seq_list_cp) == 0:
        sys.exit('empty list')
    elif len(seq_list_cp) == 1:
        merged_seq = seq_list_cp[0]
    else:
        merged_seq = seq_list_cp.pop(0)
        for each_seq in seq_list_cp:
            merged_seq = merge_overlap_seq(merged_seq, each_seq)
    return merged_seq

if __name__ == '__main__':
    mirdeep2_dir = args.mirdeep2_dir
    output_dir =args.output_dir
    novel_miRNA_file = args.novel_miRNA_file

    miRNA_seq_dict = {}
    sample_list = os.listdir(mirdeep2_dir)
    for each_sample in sample_list:
        each_sample_dir = os.path.join(mirdeep2_dir, each_sample)        
        each_novel_miRNA_seq_file = glob.glob(r'%s/result*.csv' % each_sample_dir)[0]
        with open(each_novel_miRNA_seq_file) as each_novel_miRNA_seq_file_info:
            flag = 0
            for n, eachline in enumerate(each_novel_miRNA_seq_file_info):
                if 'novel miRNAs predicted by miRDeep2' in eachline:
                    flag = 1
                if 'mature miRBase miRNAs detected by miRDeep2' in eachline:
                    flag = 0
                eachline_info = eachline.strip().split('\t')
                if flag == 1 and len(eachline_info) == 17:
                    if eachline_info[0] == 'provisional id':
                        continue
                    novel_id = ':'.join(eachline_info[16].split(':')[0:2])
                    novel_id = re.sub('\.\.','-',novel_id)
                    miRNA_seq_dict.setdefault(novel_id, {})[each_sample] = eachline_info[13:16]

    mature_fa = os.path.join(output_dir, 'novel.mature.fa')
    star_fa = os.path.join(output_dir, 'novel_star.fa')
    hairpin_fa = os.path.join(output_dir, 'novel_hairpin.fa')
    mature_fa_dict = {}
    star_fa_dict = {}
    hairpin_fa_dict = {}
    seq_file_list = [mature_fa, star_fa, hairpin_fa]
    seq_dict_list = [mature_fa_dict, star_fa_dict, hairpin_fa_dict]

    sample_index_dict = {}
    with open(novel_miRNA_file) as novel_miRNA_file_info:
        for n, eachline in enumerate(novel_miRNA_file_info):
            eachline_info = eachline.strip().split('\t')
            if n == 0:
                for m, each_col in enumerate(eachline_info):
                    if '_novel_miRNA' in each_col:
                        sample_id = each_col.split('_novel_miRNA')[0]
                        sample_index_dict[m] = sample_id                    
            else:                
                novel_id = eachline_info[0]
                novel_seq_list = [[], [], []]
                for m, each_col in enumerate(eachline_info):
                    if m in  sample_index_dict:
                        sample_id = sample_index_dict[m]
                        if each_col != '-':
                            each_col_list = each_col.split(',')
                            each_col_list = [each.strip() for each in each_col_list]
                            for each_col_info in each_col_list:
                                for i, each_seq in enumerate(miRNA_seq_dict[each_col_info][sample_id]):
                                    novel_seq_list[i].append(miRNA_seq_dict[each_col_info][sample_id][i])
                for j, each_seq_list in enumerate(novel_seq_list):
                    each_merged_seq = merge_seq_list(each_seq_list)
                    seq_dict_list[j][novel_id] = each_merged_seq

    for n, each_dict in enumerate(seq_dict_list):
        each_file = seq_file_list[n]
        python_tools.write_obj_to_file(each_dict, each_file)


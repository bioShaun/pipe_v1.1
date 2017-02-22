import sys
import os
import glob

name_map_dict = {}

with open('name_map') as name_map_inf:
    for eachline in name_map_inf:
        eachline_inf = eachline.strip().split('\t')
        name_map_dict[eachline_inf[0]] = eachline_inf[1]


file_list = os.listdir('./')

for each_file in file_list:
    if os.path.isdir(each_file):
        wgc_id = each_file.split('_')[1]
        if wgc_id not in name_map_dict:
            continue
        sample_id = name_map_dict[wgc_id]
        for n in range(1,3):
            each_md5_file = os.path.join(each_file, '{0}_combined_R{1}.fastq.gz.md5'.format(wgc_id, n))
            each_new_md5_file = os.path.join(each_file, '{0}_R{1}.fastq.gz.md5'.format(sample_id, n))
            each_fq_file = os.path.join(each_file, '{0}_combined_R{1}.fastq.gz'.format(wgc_id, n))
            each_new_fq_file = os.path.join(each_file, '{0}_R{1}.fastq.gz'.format(sample_id, n))
            md5_code = ''
            with open(each_md5_file) as each_md5_file_inf:
                for m, eachline in enumerate(each_md5_file_inf):
                    if m == 0:
                        md5_code = eachline.strip().split()[0]
            with open(each_new_md5_file, 'w') as each_new_md5_file_inf:
                each_new_md5_file_inf.write('{0}  {1}_R{2}.fastq.gz\n'.format(md5_code, sample_id, n))
            os.system('rm {0}'.format(each_md5_file))
            os.system('mv {0} {1}'.format(each_fq_file, each_new_fq_file))
        os.system('mv Sample_{0} {1}'.format(wgc_id, sample_id))

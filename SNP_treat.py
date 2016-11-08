import sys
import os

snp_file = sys.argv[1]
sample_order = sys.argv[2]
group_sample = sys.argv[3]
out_dir = sys.argv[4]

def get_group_snp_type(ref, alt, snp_type_list):
    test_snp_type_list = list(set(snp_type_list))
    out_list = []
    if len(test_snp_type_list) == 1:
        out_list = [test_snp_type_list[0]]
    else:
        if 'NA' in test_snp_type_list:
            test_snp_type_list.remove('NA')
        if len(test_snp_type_list) == 1:
            out_list = [test_snp_type_list[0]]
        else:
            test_snp_type_list2 = []
            for each in test_snp_type_list:
                test_snp_type_list2.extend(list(each))
            alt_list = alt.split(',')
            if ref in test_snp_type_list2:
                for each_alt in alt_list:
                    if each_alt in test_snp_type_list2:
                        out_list.append('{ref}{each_alt}'.format(**locals()))
            else:
                for each_alt in alt_list:
                    if each_alt in test_snp_type_list2:
                        out_list.append('{each_alt}{each_alt}'.format(**locals()))
    return out_list

def check_dgu_specific_snp(snp_list):
    dgu_snp = list(snp_list[0])[0]
    other_snp = [snp_list[1]]
    other_snp.extend(snp_list[3:])
    out_flag = True
    for each in other_snp:
        if dgu_snp in each:
            out_flag = False
            break
    return out_flag

group_snp_dict = {}
group_dict = {}
sample_group_dict = {}
with open(group_sample) as group_sample_inf:
    for eachline in group_sample_inf:
        eachline_inf = eachline.strip().split('\t')
        sample_id = eachline_inf[0]
        group_id = eachline_inf[1]
        sample_group_dict[sample_id] = group_id
        group_dict[group_id] = 1

sample_order_list = [each.strip() for each in open(sample_order)]

group_list = ['Dgu']
tmp_group_list = sorted(group_dict.keys())
group_list.extend(tmp_group_list)

all_out = os.path.join(out_dir,'all.snp.inf.txt')
dgu_out = os.path.join(out_dir,'Dgu.snp.inf.txt')
all_out_inf = open(all_out, 'w')
all_out_inf.write('Chr\tPos\tRef\tAlt\t%s\n' % '\t'.join(group_list))
dgu_out_inf = open(dgu_out, 'w')
dgu_out_inf.write('Chr\tPos\tRef\tAlt\t%s\n' % '\t'.join(group_list))

with open(snp_file) as snp_file_inf:
    for eachline in snp_file_inf:
        group_snp_dict = {}
        eachline_inf = eachline.strip().split('\t')
        each_snp_inf = eachline_inf[0:4]
        each_sample_snp_inf = eachline_inf[4:]
        ref = eachline_inf[2]
        alt = eachline_inf[3]
        haptype = '{ref}{alt}'.format(**locals())
        each_line_out = []
        for n, each_snp in enumerate(each_sample_snp_inf):
            each_snp_sample_id = sample_order_list[n]
            if 'Dgu' in each_snp_sample_id:
                each_snp_group_id = 'Dgu'
            else:
                each_snp_group_id = sample_group_dict[each_snp_sample_id]
            if each_snp_group_id not in group_snp_dict:
                group_snp_dict[each_snp_group_id] = [each_snp]
            else:
                group_snp_dict[each_snp_group_id].append(each_snp)
        for each_group in group_list:
            each_group_snp = get_group_snp_type(ref, alt, group_snp_dict[each_group])
            each_group_snp_out = ','.join(each_group_snp)
            each_line_out.append(each_group_snp_out)
        all_out_inf.write('%s\t%s\n' % ('\t'.join(each_snp_inf), '\t'.join(each_line_out)))
        if haptype not in each_line_out[0] and each_line_out[0] != "NA":
            dgu_flag = check_dgu_specific_snp(each_line_out)
            if dgu_flag:
                dgu_out_inf.write('%s\t%s\n' % ('\t'.join(each_snp_inf), '\t'.join(each_line_out)))
all_out_inf.close()
dgu_out_inf.close()

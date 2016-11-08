import sys
import os
import numpy
import json
import python_tools

if not len(sys.argv) == 6:
    print '    python ' + sys.argv[0] + ' group_sample exp.table out.dir exp.cutoff rep.cutoff'
    sys.exit(0)

group_sample = sys.argv[1]
exp_table = sys.argv[2]
out_dir = sys.argv[3]

group_dict = {}
sample_dict = {}
group_exp_dict = {}
sample_index_dict = {}
reproducibility_dict = {}
tpm_cutoff = (0.5, 1, 2, 4, 8)

EXP_CUTOFF = float(sys.argv[4])
REP_CUTOFF = float(sys.argv[5])

def write_list_to_file(obj, fn, header = ''):
    fh = open(fn, 'w')
    if header != '':
        fh.write('%s\n' % header)
    for each_list in obj:
        each_list = [str(each) for each in each_list]
        eachline = '\t'.join(each_list)
        fh.write('%s\n' % eachline)
    fh.close() 

def get_group_reproducibility(gene, group_exp_dict, group_dict, reproducibility_dict): 
    non_rep_sample_list = []
    non_rep_group_list = [] 
    tpm_max_list = []
    for each_group in group_exp_dict[gene]:
        sample_list = group_dict[each_group][0]
        exp_list = group_exp_dict[gene][each_group]
        tpm_max_list.append(max(exp_list))              
        exp_status_list = []
        for each_exp in exp_list:
            if each_exp <= EXP_CUTOFF:
                exp_status_list.append(0)
            else:
                exp_status_list.append(1)
        exp_status = sum(exp_status_list)
#        if exp_status == 0 or exp_status == len(exp_status_list):
#        ref_flag = 0
        if len(exp_status_list) <= 3:
            if exp_status == 0 or exp_status == len(exp_status_list):
                pass
            else:
                non_rep_group_list.append(each_group)
        else:
            if (exp_status > len(exp_status_list)*REP_CUTOFF) or (exp_status < len(exp_status_list)*(1-REP_CUTOFF)):
                pass
            else:
                non_rep_group_list.append(each_group)
        if exp_status < (len(exp_status_list)/2.0):
            for n, each_stat in enumerate(exp_status_list):
                if each_stat == 1:
                    try:
                        non_rep_sample_list.append(sample_list[n])
                    except:
                        print each_group
        elif exp_status > (len(exp_status_list)/2.0):
            for n, each_stat in enumerate(exp_status_list):
                if each_stat == 0:
                    try:
                        non_rep_sample_list.append(sample_list[n])
                    except:
                        print each_group
        else:
            non_rep_sample_list.append('NS')
    tpm_max = max(tpm_max_list)
    for m, each_cut in enumerate(tpm_cutoff):
        if tpm_max < each_cut:
            if m == 0:
                break
            else:
                reproducibility_dict.setdefault(m-1, {})[gene] = [non_rep_group_list, non_rep_sample_list]
                break
    else:
        reproducibility_dict.setdefault(m, {})[gene] = [non_rep_group_list, non_rep_sample_list]
    return non_rep_group_list, non_rep_sample_list, tpm_max

## 
with open(group_sample) as group_sample_inf:
    for eachline in group_sample_inf:
        eachline_inf = eachline.strip().split('\t')
        group_name = eachline_inf[0]
        sample_name = eachline_inf[1]
        if group_name not in group_dict:
            group_dict[group_name] = [[sample_name], []]
        else:
            group_dict[group_name][0].append(sample_name)
        sample_dict[sample_name] = group_name

group_num = len(group_dict.keys())
all_non_rep_group_list = []
all_non_rep_sample_list = []

## out put file names
gene_rep_status_file = os.path.join(out_dir, 'gene.reproducibility.txt')
gene_rep_detail_file = os.path.join(out_dir, 'gene.exp.reproducibility.txt')
gene_rep_summary_file = os.path.join(out_dir, 'gene.exp.summary.txt')
sample_file = os.path.join(out_dir, 'sample_reproducibility.txt')
group_file = os.path.join(out_dir, 'group_reproducibility.txt')
group_file_plot = os.path.join(out_dir, 'group_reproducibility_plot.txt')
gene_variability_plot_file = os.path.join(out_dir, 'gene.var.plot.txt')
gene_variability_file = os.path.join(out_dir, 'gene.variability.txt')
group_exp_dict_json = os.path.join(out_dir, 'gene.exp.dict.json')
gene_rep_status = open(gene_rep_status_file, 'w')
gene_rep_status.write('Gene\treproducibility_stat\treproducibility_percentage\tnon_reproducible_group\tnon_reproducible_sample\n')
all_exp_list = []

with open(exp_table) as exp_table_inf:
    for n, eachline in enumerate(exp_table_inf):
        eachline_inf = eachline.strip().split('\t')
        if n == 0:
            for m, each_sample in enumerate(eachline_inf):
                if m > 0:
                    sample_index_dict[m] = each_sample
        else:
            gene_id = eachline_inf[0]
            for m, each_exp in enumerate(eachline_inf):
                eachline_exp_inf = [float(each) for each in eachline_inf[1:]]
                all_exp_list.append(numpy.mean(eachline_exp_inf))
                if m > 0:
                    sample_name = sample_index_dict[m]
                    group_name = sample_dict[sample_name]
                    sample_num = len(group_dict[group_name][0])
                    group_dict_index = group_dict[group_name][0].index(sample_name)                    
                    if gene_id in group_exp_dict and group_name in group_exp_dict[gene_id]:
                        group_exp_dict[gene_id][group_name][group_dict_index] = float(eachline_inf[m])
                    else:
                        group_exp_dict.setdefault(gene_id, {})[group_name] = [''] * (sample_num)
                        group_exp_dict[gene_id][group_name][group_dict_index] = float(eachline_inf[m])                
            for each_group in group_exp_dict[gene_id]:
                each_group_exp_list = group_exp_dict[gene_id][each_group][:]
                each_group_exp_list.insert(0, gene_id)
                group_dict[each_group][1].append(each_group_exp_list)
            non_rep_group_list, non_rep_sample_list, tpm_max = get_group_reproducibility(gene_id, group_exp_dict, group_dict, reproducibility_dict)
            all_non_rep_group_list.extend(non_rep_group_list)
            all_non_rep_sample_list.extend(non_rep_sample_list)
            rep_num = group_num - len(non_rep_group_list)
            rep_stat = '%s/%s' % (rep_num, group_num)
            rep_percentage = round(100*rep_num/float(group_num), 2)
            non_rep_group_out = ','.join(non_rep_group_list)
            non_rep_sample_out = ','.join(non_rep_sample_list)
            gene_rep_status.write('{gene_id}\t{rep_stat}\t{rep_percentage}\t{non_rep_group_out}\t{non_rep_sample_out}\t{tpm_max}\n'.format(**locals()))

if not os.path.exists(group_exp_dict_json):
    python_tools.write_obj_to_json(group_exp_dict, group_exp_dict_json)

all_exp_mean = numpy.mean(all_exp_list)

gene_count = len(group_exp_dict.keys())
gene_rep_summary_list = []
tpm_breaks = ['Group']
tmp_breaks_genes = []
group_rep_dict = {}
with open(gene_rep_detail_file, 'w') as gene_rep_detail:
    for n, each_tpm in enumerate(tpm_cutoff):
        non_rep_group = []
        non_rep_sample = []        
        if n+1 < len(tpm_cutoff):
            flag = '%s<TPM<=%s' % (tpm_cutoff[n], tpm_cutoff[n+1])
        else:
            flag = 'TPM>%s' % tpm_cutoff[n]
        tpm_breaks.append(flag)
        gene_rep_summary_list.append([flag, 0, 0, 0])
        for each_gene in reproducibility_dict[n]:
            each_gene_non_rep_list = reproducibility_dict[n][each_gene][0]
            non_rep_group_num = len(reproducibility_dict[n][each_gene][0])
            non_rep_group.extend(reproducibility_dict[n][each_gene][0])
            non_rep_sample.extend(reproducibility_dict[n][each_gene][1])
            rep_percentage = 1-(non_rep_group_num/float(group_num))
            gene_rep_summary_list[n][1] += 1
            if rep_percentage >= REP_CUTOFF:
                rep_stat = 1
                gene_rep_summary_list[n][2] += 1
            else:
                rep_stat = 0  
            rep_percentage_out = round(rep_percentage*100, 2)
            gene_rep_detail.write('%s\t%s\t%s\t%s\n' % (each_gene, rep_stat, rep_percentage_out, flag))
            for each_group in each_gene_non_rep_list:
                if each_group in group_rep_dict:
                    group_rep_dict[each_group][n] += 1
                else:
                    group_rep_dict[each_group] = [0] * len(tpm_cutoff)
                    group_rep_dict[each_group][n] += 1
        each_tpm_rep_percentage = round(gene_rep_summary_list[n][2]/float(gene_rep_summary_list[n][1]), 2)
        tmp_breaks_genes.append(gene_rep_summary_list[n][1])
        gene_rep_summary_list[n][3] = each_tpm_rep_percentage

write_list_to_file(gene_rep_summary_list, gene_rep_summary_file)

sample_file_inf = open(sample_file, 'w')
group_file_inf = open(group_file, 'w')
group_file_plot_inf = open(group_file_plot, 'w')
group_file_inf.write('%s\n' % '\t'.join(tpm_breaks))
group_file_plot_inf.write('Group\texp\tcount\tprecent\n')
for each_group in group_dict:
    for each_sample in group_dict[each_group][0]:
        each_sample_non_rep_count = non_rep_sample.count(each_sample)
        each_sample_non_rep_percent = round(each_sample_non_rep_count/float(gene_count), 2)
        sample_file_inf.write('%s\t%s\t%s\t%s\n' % (each_group, each_sample, each_sample_non_rep_count, each_sample_non_rep_percent))
    if each_group not in group_rep_dict:
        group_rep_dict[each_group] = [0] * len(tpm_cutoff)
    each_group_rep_list = []
    for n, each_no_rep_num  in enumerate(group_rep_dict[each_group]):
        each_tpm_genes = tmp_breaks_genes[n]
        each_group_rep_gens = each_tpm_genes - each_no_rep_num
        each_group_rep_rate = '%s%%' % round(100*each_group_rep_gens/float(each_tpm_genes), 2)
        each_group_rep_rate_num = round(each_group_rep_gens/float(each_tpm_genes), 2)
        each_group_rep_list.append('%s(%s)' % (each_group_rep_gens, each_group_rep_rate))
        group_file_plot_inf.write('%s\t%s\t%s\t%s\n' % (each_group, tpm_breaks[n+1], each_group_rep_gens, each_group_rep_rate_num ))
    group_file_inf.write('%s\t%s\n' % (each_group, '\t'.join(each_group_rep_list)))
group_file_plot_inf.close()
group_file_inf.close()    
sample_file_inf.close()

## gene variability
group_list = group_dict.keys()
gene_variability_file_inf = open(gene_variability_file, 'w')
gene_variability_file_inf.write('Gene\t%s\n' % '\t'.join(group_list))
gene_variability_plot_file_inf = open(gene_variability_plot_file, 'w')
gene_variability_plot_file_inf.write('gene\tgroup\tvar\n')
for each_gene in group_exp_dict:
    each_gene_var_list = []
    for each_group in group_list:
        each_group_exp = group_exp_dict[each_gene][each_group]
        if max(each_group_exp) > 0.2:
            each_var = str(numpy.std(each_group_exp)/all_exp_mean)
            gene_variability_plot_file_inf.write('%s\t%s\t%s\n' % (each_gene, each_group, each_var))
        else:
            each_var = 'NA'
        each_gene_var_list.append(each_var)
    gene_variability_file_inf.write('%s\t%s\n' % (each_gene, '\t'.join(each_gene_var_list)))
gene_variability_file_inf.close()
gene_variability_plot_file_inf.close()

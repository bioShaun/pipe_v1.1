'''
get kegg pathview plots

'''

import sys
import os
import argparse
import python_tools

parser = argparse.ArgumentParser()
parser.add_argument('--kegg_table', help = 'KOBAS kegg enrichment analysis result', required = True)
parser.add_argument('--blast_out', help = 'ko blast result', required = True)
parser.add_argument('--species', help = 'kegg species', required = True)
parser.add_argument('--diff_out', help = 'diff analysis result', required = True)
parser.add_argument('--diff_method', help = 'diff analysis method',choices = ['edgeR','DESeq2','Cuffdiff'], default = 'edgeR')
parser.add_argument('--out_dir', help = 'output directory', required = True)
args = parser.parse_args()

pathview_script = '/home/lxgui/scripts/kegg_pathview_plot.R'
kegg_pathway_dir = '/home/public/database/kegg/pathway/'

exclude_pathway_id = ('01100','01110')
problem_pathway_id= ('bta04215')

def get_fc_row(diff_method):
    if diff_method == 'edgeR' :
        return 2
    elif diff_method == 'Cuffdiff' :
        return 6
    elif diff_method == 'DESeq2' :
        return 5
    else :
        sys.exit('wrong diff analysis method! : %s ' % diff_method)

def get_diff_fc_dict(diff_out,diff_method) :
    fc_row = get_fc_row(diff_method)
    diff_fc_dict = {}
    # fc_list = []
    # inf_flag = 0
    with open(diff_out) as diff_out_info :
        for n,eachline in enumerate(diff_out_info) :
            if n != 0 :
                eachline_info = eachline.strip().split('\t')
                if eachline_info[fc_row-1] != 'NA' :
                    qurey_id = eachline_info[0]
                    fc = float(eachline_info[fc_row-1])
                    diff_fc_dict[qurey_id] = fc
    return diff_fc_dict

def get_kegg_map(kegg_blast_out) :
    kegg_map_dict = {}
    kegg_to_gene_dict = {}
    kegg_identity_dict = {}
    with open(kegg_blast_out) as kegg_blast_out_info :
        for eachline in kegg_blast_out_info :
            eachline_info = eachline.strip().split('\t')
            gene_id = eachline_info[0]
            kegg_gene = eachline_info[1].split(':')[1]
            identity = float(eachline_info[3])
            if kegg_gene not in kegg_to_gene_dict :
                kegg_to_gene_dict[kegg_gene] = gene_id
                kegg_identity_dict[kegg_gene] = identity
            else :
                if identity > kegg_identity_dict[kegg_gene] :
                    kegg_to_gene_dict[kegg_gene] = gene_id
                    kegg_identity_dict[kegg_gene] = identity
    kegg_map_dict = python_tools.invert_dict(kegg_to_gene_dict)
    return kegg_map_dict

def plot_pathview(species,pathview_id,each_pathway_kegg_fc_out,out_dir):
    cmd = 'Rscript %s %s %s %s %s' % (pathview_script,species,pathview_id,each_pathway_kegg_fc_out,out_dir)
    python_tools.circ_call_process(cmd)
    os.system('rm %s' % each_pathway_kegg_fc_out)

def kegg_pathway_plot(kegg_out,kegg_map_dict,diff_fc_dict):
    with open(kegg_out) as kegg_out_info :
        for eachline in kegg_out_info :            
            eachline_info = eachline.strip().split('\t')
            if (not eachline_info[0].startswith('#')) and len(eachline_info) == 9 :
                pathway_id = eachline_info[2]
                pathview_id = pathway_id.split(args.species)[1]
                kegg_pic1 = os.path.join(args.out_dir,'%s.pathview.png' % pathway_id)
                kegg_pic2 = os.path.join(args.out_dir,'%s.pathview.gene.png' % pathway_id)
                if os.path.isfile(kegg_pic1) and os.stat(kegg_pic1).st_size and os.path.isfile(kegg_pic2) and os.stat(kegg_pic2).st_size:
                    pass
                elif pathview_id in exclude_pathway_id or pathway_id in problem_pathway_id:
                    target_pic = os.path.join(kegg_pathway_dir,'%s/%s.png' % (args.species, pathway_id))
                    os.system('cp %s %s' % (target_pic,kegg_pic1))
                    os.system('cp %s %s' % (target_pic,kegg_pic2)) 
                else :
                    each_pathway_kegg_fc = os.path.join(args.out_dir,'%s.keggid' % pathway_id)
                    each_pathway_kegg_fc_out = open(each_pathway_kegg_fc,'w')
                    pathway_gene_list = eachline_info[7].split('|')
                    for each_gene in pathway_gene_list :
                        if each_gene in kegg_map_dict :
                            kegg_id = kegg_map_dict[each_gene]
                            fc = diff_fc_dict[each_gene]
                            each_pathway_kegg_fc_out.write('%s\t%s\n' % (kegg_id,fc))
                    each_pathway_kegg_fc_out.close()
                    plot_pathview(args.species,pathview_id,each_pathway_kegg_fc,args.out_dir)

if __name__ == '__main__':
    diff_fc_dict = get_diff_fc_dict(args.diff_out,args.diff_method)
    kegg_map_dict = get_kegg_map(args.blast_out)
    kegg_pathway_plot(args.kegg_table,kegg_map_dict,diff_fc_dict)



import sys
import RNAseq_tools
from scipy.stats import hypergeom
import os

if not len(sys.argv) == 7:
    print '    python ' + sys.argv[0] + ' motif_out network_edge_file network_node_file gtf out_dir promoter_sum'
    sys.exit(0)

motif_out = sys.argv[1]
network_edge_file = sys.argv[2]
network_node_file = sys.argv[3]
gtf = sys.argv[4]
out_dir = sys.argv[5]
promoter_sum = int(sys.argv[6])

motif_dict = {}
gene_dict = {}
tr_dict = RNAseq_tools.get_transcript_info(gtf)

def get_hypergeom_p(k, M, n, N):
    rv = hypergeom(M, n, N)
    pvalue = 0
    for m in range(k + 1, n + 1):
        pvalue += rv.pmf(m)
    return pvalue

with open(motif_out) as motif_out_info:
    for n, eachline in enumerate(motif_out_info):
        if n != 0:
            eachline_info = eachline.strip().split('\t')
            tr_id = eachline_info[1]
            gene_id = tr_dict[tr_id]['gene_id']
            motif_id = eachline_info[0]
            dist = int(eachline_info[2]) - 1000
            motif_len = len(eachline_info[-1])        
            if motif_id not in motif_dict:
                motif_dict.setdefault(motif_id, {})['gene'] = set([gene_id])
                motif_dict.setdefault(motif_id, {})['len'] = motif_len
            else:
                motif_dict[motif_id]['gene'].add(gene_id)
            if gene_id in gene_dict:                                
                if motif_id in gene_dict[gene_id]['motif_id']:
                    motif_index = gene_dict[gene_id]['motif_id'].index(motif_id)
                    if dist > gene_dict[gene_id]['motif_dis'][motif_index]:
                        gene_dict[gene_id]['motif_dis'][motif_index] = dist
                else:
                    gene_dict.setdefault(gene_id, {})['motif_id'].append(motif_id)
                    gene_dict.setdefault(gene_id, {})['motif_dis'].append(dist)
            else:
                gene_dict.setdefault(gene_id, {})['motif_id'] = [motif_id]
                gene_dict.setdefault(gene_id, {})['motif_dis'] = [dist]

with open(network_edge_file) as network_edge_file_info:
    for n, eachline in enumerate(network_edge_file_info):
        if n != 0:
            eachline_info = eachline.strip().split('\t')
            fromNode = eachline_info[0]
            toNode = eachline_info[1]
            if fromNode in gene_dict and 'net' in gene_dict[fromNode]:
                gene_dict[fromNode]['net'].append(toNode)
            else:
                gene_dict.setdefault(fromNode, {})['net'] = [toNode]
                
with open(network_node_file) as network_node_file_info:
    for n, eachline in enumerate(network_node_file_info):
        if n != 0:
            eachline_info = eachline.strip().split('\t')
            nodeName = eachline_info[0]
            nodecol = eachline_info[-1]
            gene_dict.setdefault(nodeName, {})['col'] = [nodecol]

motif_num_enrich_file = os.path.join(out_dir, 'motif_enrich2.txt')

with open(motif_num_enrich_file, 'w') as motif_num_enrich_file_info:
    for each_gene in gene_dict:
        if 'motif_id' not in gene_dict[each_gene]:
            continue
        for m, each_motif in enumerate(gene_dict[each_gene]['motif_id']):
            each_motif_total = len(motif_dict[each_motif]['gene'])
            each_motif_len = motif_dict[each_motif]['len']
            each_motif_num = 1
            each_motif_dist = [gene_dict[each_gene]['motif_dis'][m]]
            each_promoter_num = 1
            if 'MYB' not in each_motif and 'WRKY' not in each_motif:
                continue
            if 'net' not in gene_dict[each_gene]:
                continue
            for each_con in gene_dict[each_gene]['net']:
                each_promoter_num += 1
                if each_con not in gene_dict:
                    continue
                if 'motif_id' not in gene_dict[each_con]:
                    continue
                if each_motif in gene_dict[each_con]['motif_id']:
                    each_motif_num += 1
                    each_motif_index = gene_dict[each_con]['motif_id'].index(each_motif)
                    each_motif_dist.append(gene_dict[each_con]['motif_dis'][each_motif_index])
            average_dis = sum(each_motif_dist)/float(len(each_motif_dist))
            tmp1 = ((1000 - each_motif_len + 1)**2 - 1)/float(len(each_motif_dist))
            zscore = (500 + average_dis)/(tmp1**0.5)
            pvalue = get_hypergeom_p(each_motif_num, promoter_sum, each_motif_total, each_promoter_num)
            motif_num_enrich_file_info.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (each_gene, each_motif, each_motif_total, promoter_sum, each_promoter_num, each_motif_num, pvalue, zscore))








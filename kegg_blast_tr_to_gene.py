import sys
import os
import RNAseq_tools
import python_tools

if not len(sys.argv) == 4:
    print '    python ' + sys.argv[0] + ' gtf.file tr.blast.annotate \
gene.blast.annotate'
    sys.exit(0)

gtf_file = sys.argv[1]
tr_blast = sys.argv[2]
gene_blast = sys.argv[3]

tr_dict = RNAseq_tools.get_transcript_info(gtf_file)

gene_dict = {}
ko_dict = {}
out_dict = {}
out_inf_list = []

with open(tr_blast) as tr_blast_inf:
    for eachline in tr_blast_inf:
        eachline_inf = eachline.strip().split('\t')
        gene_id = tr_dict[eachline_inf[0]]['gene_id']
        ko_id = eachline_inf[1]
        bitscore = float(eachline_inf[-1])
        eachline_out = '%s\t%s\n' % (gene_id, '\t'.join(eachline_inf[1:]))
        if gene_id not in gene_dict and ko_id not in ko_dict:
            gene_dict[gene_id] = [ko_id, bitscore]
            ko_dict[ko_id] = [gene_id, bitscore]
            out_dict.setdefault(gene_id, {})[ko_id] = eachline_out
        elif gene_id in gene_dict:
            if bitscore > gene_dict[gene_id][-1]:
                old_ko = gene_dict[gene_id][0]
                del ko_dict[old_ko]
                del gene_dict[gene_id]
                del out_dict[gene_id][old_ko]
                if ko_id in ko_dict:
                    if bitscore > ko_dict[ko_id][-1]:
                        old_gene = ko_dict[ko_id][0]
                        del ko_dict[ko_id]
                        del gene_dict[old_gene]
                        del out_dict[old_gene][ko_id]
                        ko_dict[ko_id] = [gene_id, bitscore]
                        gene_dict[gene_id] = [ko_id, bitscore]
                        out_dict.setdefault(gene_id, {})[ko_id] = eachline_out
                else:
                    ko_dict[ko_id] = [gene_id, bitscore]
                    gene_dict[gene_id] = [ko_id, bitscore]
                    out_dict.setdefault(gene_id, {})[ko_id] = eachline_out
        else:
            if bitscore > ko_dict[ko_id][-1]:
                old_gene = ko_dict[ko_id][0]
                del out_dict[old_gene][ko_id]
                del gene_dict[old_gene]
                ko_dict[ko_id] = [gene_id, bitscore]
                gene_dict[gene_id] = [ko_id, bitscore]
                out_dict.setdefault(gene_id, {})[ko_id] = eachline_out

with open(gene_blast, 'w') as gene_blast_inf:
    for each_gene in out_dict:
        for each_ko in out_dict[each_gene]:
            gene_blast_inf.write(out_dict[each_gene][each_ko])

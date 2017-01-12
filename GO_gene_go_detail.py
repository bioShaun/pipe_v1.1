import sys
import os

if not len(sys.argv) == 4:
    print '    python ' + sys.argv[0] + ' gene.go go.detail gene.go.detail'
    sys.exit(0)

gene_go = sys.argv[1]
go_detail = sys.argv[2]
gene_go_detail = sys.argv[3]

gene_go_dict = {}
go_inf_dict = {}

with open(go_detail) as go_detail_inf:
    for n, eachline in enumerate(go_detail_inf):
        if n != 0:
            eachline_inf = eachline.strip().split('\t')
            go_id = eachline_inf[0]
            go_inf = eachline_inf[1:]
            if go_inf[1] != "NA":
                go_inf_dict[go_id] = go_inf

with open(gene_go) as gene_go_inf:
    for n, eachline in enumerate(gene_go_inf):
        if n != 0:
            eachline_inf = eachline.strip().split(',')
            gene_id = eachline_inf[0]
            go_id = eachline_inf[1]
            if go_id == "":
                continue
            if go_id not in go_inf_dict:
                continue
            go_ontology = go_inf_dict[go_id][0]
            go_term = go_inf_dict[go_id][1]
            if gene_id not in gene_go_dict:
                gene_go_dict[gene_id] = [[go_id], [go_ontology], [go_term]]
            else:
                gene_go_dict[gene_id][0].append(go_id)
                gene_go_dict[gene_id][1].append(go_ontology)
                gene_go_dict[gene_id][2].append(go_term)

with open(gene_go_detail, 'w') as gene_go_detail_inf:
    gene_go_detail_inf.write('Gene_ID\tGO_ID\tGO_Ontology\tGO_Term\n')
    for each_gene in gene_go_dict:
        go_id = '||'.join(gene_go_dict[each_gene][0])
        go_ontology = '||'.join(gene_go_dict[each_gene][1])
        go_term = '||'.join(gene_go_dict[each_gene][2])
        gene_go_detail_inf.write('{each_gene}\t{go_id}\t{go_ontology}\t{go_term}\n'.format(**locals()))


import sys
import os
import python_tools

if not len(sys.argv) == 4:
    print '    python ' + sys.argv[0] + ' gene.interpro.id interpro.map gene.interpro.des'
    sys.exit(0)

gene_interpro_id_file = sys.argv[1]
interpro_map = sys.argv[2]
gene_interpro_des_file = sys.argv[3]

interpro_map_dict = {}

if interpro_map.endswith('json'):
    interpro_map_dict = python_tools.load_fn_to_obj(interpro_map)
else:
    with open(interpro_map) as interpro_map_inf:
        for eachline in interpro_map_inf:
            eachline_inf = eachline.strip().split('\t')
            interpro_map_dict[eachline_inf[0]] = eachline_inf[1]
    interpro_map_json = '%s.json' % interpro_map
    python_tools.write_obj_to_json(interpro_map_dict, interpro_map_json)

gene_inf_dict = {}
with open(gene_interpro_id_file) as gene_interpro_id:
    for n, eachline in enumerate(gene_interpro_id):
        if n != 0:
            eachline_inf = eachline.strip().split(',')
            gene_id = eachline_inf[0]
            gene_name = eachline_inf[1]
            interpro_id = eachline_inf[2]            
            if gene_id not in gene_inf_dict:
                gene_inf_dict[gene_id] = [[], [], []]
            if gene_name != '':
                if gene_name not in gene_inf_dict[gene_id][0]:
                    gene_inf_dict[gene_id][0].append(gene_name)
            if interpro_id != '':
                if interpro_id not in interpro_map_dict:
                    continue
                interpro_des = interpro_map_dict[interpro_id]
                gene_inf_dict[gene_id][1].append(interpro_id)
                gene_inf_dict[gene_id][2].append(interpro_des)
python_tools.write_obj_to_json(gene_inf_dict, gene_interpro_des_file)


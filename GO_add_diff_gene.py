import os
import sys

if not len(sys.argv) == 5:
    print '    python ' + sys.argv[0] + ' go.table diff.id go.id go.table.annotate'
    sys.exit(0)

go_table = sys.argv[1]
diff_id = sys.argv[2]
go_id = sys.argv[3]
go_anno = sys.argv[4]

diff_id_list = [each.strip() for each in open(diff_id)]

go_dict = {}
with open(go_id) as go_id_inf:
    for n, eachline in enumerate(go_id_inf):
        if n != 0:
            eachline_inf = eachline.strip().split(',')
            gene_id = eachline_inf[0]
            go_id = eachline_inf[1]
            if gene_id in diff_id_list and go_id != '':
                if go_id in go_dict :
                    if gene_id not in go_dict[go_id]:
                        go_dict[go_id].append(gene_id)
                else:
                    go_dict[go_id] = [gene_id]

go_anno_inf = open(go_anno, 'w')
with open(go_table) as go_table_inf:
    for n, eachline in enumerate(go_table_inf):
        eachline = eachline.strip()
        if n == 0:
            go_anno_inf.write('%s\tDE_id\n' % eachline)
        else:
            eachline_inf = eachline.split('\t')
            go_id = eachline_inf[0]
            de_genes = ','.join(go_dict[go_id])
            numDEInCat = int(eachline_inf[3])
            numDEInCat2 = len(go_dict[go_id])
            if numDEInCat != numDEInCat2:
                sys.exit('numDEInCat in %s not equal in topGO(%s) and this program(%s)!' % (go_id, numDEInCat, numDEInCat2))
            go_anno_inf.write('%s\t%s\n' % (eachline, de_genes))
go_anno_inf.close()

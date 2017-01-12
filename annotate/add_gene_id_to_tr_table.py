import pandas as pd
import sys
import os

if not len(sys.argv) == 5:
    print '    python ' + sys.argv[0] + ' tr.table gene.tr.map tr.column tr.gene.table'
    sys.exit(1)

tr_table_path = sys.argv[1]
gene_tr_map = sys.argv[2]
tr_column = int(sys.argv[3])
tr_gene_table_path = sys.argv[4]

tr_table_df = pd.read_table(tr_table_path, index_col = 0)
gene_tr_map_df = pd.read_table(gene_tr_map, header = None, index_col = 1)
gene_tr_map_df.columns = ["Gene_ID"]
if tr_column == 1:
    add_gene_table = pd.merge(tr_table_df, gene_tr_map_df, left_index = True, right_index = True)
else:
    tr_column_name = tr_table_df.columns[(tr_column-2)]
    add_gene_table = pd.merge(tr_table_df, gene_tr_map_df, left_on = tr_column_name, right_index = True)
sort_add_gene_table = add_gene_table.sort_index()
sort_add_gene_table.to_csv(tr_gene_table_path, sep = "\t")

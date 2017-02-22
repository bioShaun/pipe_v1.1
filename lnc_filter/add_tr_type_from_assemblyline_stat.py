import pandas as pd
import sys

input_data = sys.argv[1]
nc_list_file = sys.argv[2]
out_data = sys.argv[3]


input_df = pd.read_table(input_data, sep = '\t', index_col = 0)
nc_list = [each.strip() for each in open(nc_list_file)]

GENCODE_CATEGORY_MAP = {'IG_C_gene': 'protein_coding',
                     'IG_D_gene': 'protein_coding',
                     'IG_J_gene': 'protein_coding',
                     'IG_V_gene': 'protein_coding',
                     'TR_C_gene': 'protein_coding',
                     'TR_J_gene': 'protein_coding',
                     'TR_V_gene': 'protein_coding',
                     'TR_D_gene': 'protein_coding',
                     'TEC': 'protein_coding',
                     'nonsense_mediated_decay': 'protein_coding',
                     'non_stop_decay': 'protein_coding',
                     'retained_intron': 'protein_coding',
                     'protein_coding': 'protein_coding',
                     'ambiguous_orf': 'protein_coding',
                     'processed_transcript': 'protein_coding',
                     'Mt_rRNA': 'lncRNA',
                     'Mt_tRNA': 'lncRNA',
                     'miRNA': 'lncRNA',
                     'misc_RNA': 'lncRNA',
                     'rRNA': 'lncRNA',
                     'snRNA': 'lncRNA',
                     'snoRNA': 'lncRNA',
                     '3prime_overlapping_ncrna': 'lncRNA',
                     'lincRNA': 'lncRNA',
                     'sense_intronic': 'lncRNA',
                     'sense_overlapping': 'lncRNA',
                     'antisense': 'lncRNA',
                     'IG_C_pseudogene': 'pseudogene',
                     'IG_J_pseudogene': 'pseudogene',
                     'IG_V_pseudogene': 'pseudogene',
                     'TR_V_pseudogene': 'pseudogene',
                     'TR_J_pseudogene': 'pseudogene',
                     'Mt_tRNA_pseudogene': 'pseudogene',
                     'tRNA_pseudogene': 'pseudogene',
                     'snoRNA_pseudogene': 'pseudogene',
                     'snRNA_pseudogene': 'pseudogene',
                     'scRNA_pseudogene': 'pseudogene',
                     'rRNA_pseudogene': 'pseudogene',
                     'misc_RNA_pseudogene': 'pseudogene',
                     'miRNA_pseudogene': 'pseudogene',
                     'pseudogene': 'pseudogene',
                     'processed_pseudogene': 'pseudogene',
                     'polymorphic_pseudogene': 'pseudogene',
                     'retrotransposed': 'pseudogene',
                     'transcribed_processed_pseudogene': 'pseudogene',
                     'transcribed_unprocessed_pseudogene': 'pseudogene',
                     'unitary_pseudogene': 'pseudogene',
                     'unprocessed_pseudogene': 'pseudogene'}


tr_type_list = []
for each_rowname in input_df.index:
    if input_df.loc[each_rowname]['category'] == 'read_through':
        tr_type_list.append('read_through')
    elif input_df.loc[each_rowname]['category'] == 'same_strand':
        ref_gene_type = input_df.loc[each_rowname]['ref_gene_type']
        tr_type_list.append(GENCODE_CATEGORY_MAP[ref_gene_type])
    else:
        if each_rowname in nc_list:
            tr_type_list.append('lncRNA')
        else:
            tr_type_list.append('TUCP')

input_df['transcript_type'] = tr_type_list

input_df.to_csv(out_data, sep = '\t')

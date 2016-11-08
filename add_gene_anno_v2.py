import sys
import os

if not len(sys.argv) == 4:
    print '    python ' + sys.argv[0] + ' gene.table gene.anno.files(sep with ",") gene.table.anno'
    sys.exit(1)

gene_table = sys.argv[1]
gene_anno_file_list = sys.argv[2].strip().split(',')
gene_table_anno = sys.argv[3]

def annotate_table(table, anno, output):
    '''
    annotate gene description to the end of a table
    '''
    header = ""
    gene_anno_dict = {}
    anno_col = 0

    with open(anno) as gene_anno_inf:
        for n, eachline in enumerate(gene_anno_inf):
            eachline_inf = eachline.strip().split('\t')
            if n == 0:
                anno_col = len(eachline_inf[1:])
                header = '\t'.join(eachline_inf[1:])
            else:
                gene_id = eachline_inf[0]
                anno = '\t'.join(eachline_inf[1:])
                gene_anno_dict[gene_id] = anno

    gene_table_anno_inf = open(output, 'w')
    with open(table) as gene_table_inf:
        for n, eachline in enumerate(gene_table_inf):
            eachline_inf = eachline.strip().split('\t')
            gene_id = eachline_inf[0]
            gene_inf = '\t'.join(eachline_inf[1:])
            if n == 0:
                anno_inf = header
            else:
                if gene_id in gene_anno_dict:
                    anno_inf = gene_anno_dict[gene_id]
                else:
                    anno_inf_list = ['--'] * anno_col
                    anno_inf = '\t'.join(anno_inf_list)
            gene_table_anno_inf.write('{gene_id}\t{gene_inf}\t{anno_inf}\n'.format(**locals()))
    gene_table_anno_inf.close()

if __name__ == '__main__':
    for n, each in enumerate(gene_anno_file_list):
        each_anno_out = "{gene_table_anno}.tmp.{n}".format(**locals())
        if n+1 == len(gene_anno_file_list):
            each_anno_out = gene_table_anno
        m = n-1
        if n == 0:
            each_table = gene_table
        else:
            each_table = "{gene_table_anno}.tmp.{m}".format(**locals())
        annotate_table(each_table, each, each_anno_out)
        if n != 0:
            os.system("rm {each_table}".format(**locals()))


import sys
import os
import RNAseq_tools
import python_tools

if not len(sys.argv) == 3:
    print '    python ' + sys.argv[0] + ' gtf.file gene.length'
    sys.exit(0)

gtf_file = sys.argv[1]
gene_length = sys.argv[2]


def get_target_length_table(gtf_file, gene_length):
    if os.path.isfile(gene_length) and os.stat(gene_length).st_size :
        pass
    else :
        target_length_info = open(gene_length,'w')
        transcript_dict = RNAseq_tools.get_transcript_info(gtf_file)
        gene_dict = RNAseq_tools.get_gene_info(transcript_dict)
        for each_gene in gene_dict :
            each_gene_len = python_tools.Median(gene_dict[each_gene]['transcript_len'])
            target_length_info.write('%s\t%s\n' % (each_gene, each_gene_len))
        target_length_info.close()

if __name__ == '__main__':
    get_target_length_table(gtf_file, gene_length)
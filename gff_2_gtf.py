from HTSeq import GFF_Reader
import sys
import os

if not len(sys.argv) == 2:
    print 'python ' + sys.argv[0] + ' gff_file'
    sys.exit(0)


gff = sys.argv[1]

gff_dir, gff_name = os.path.split(gff)
gff_name = os.path.splitext(gff_name)[0]
gtf_file = os.path.join(gff_dir,'%s.gtf' % gff_name)


tr_gene_dict = {}
gtf_file_info = open(gtf_file,'w')

for eachline in GFF_Reader(gff):
    chrom = eachline.iv.chrom
    source = eachline.source
    start = eachline.iv.start + 1
    end = eachline.iv.end
    strand = eachline.iv.strand
    pos_info = '\t'.join([chrom, source, eachline.type, str(start), str(end) ,'.', strand, '.'])
    if eachline.type == 'gene':
        gene_id = eachline.attr['Name']
        tr_id = eachline.attr['Name']            
    elif eachline.type == 'mRNA':
        gene_id = eachline.attr['Parent']
        tr_id = eachline.attr['Name']
        tr_gene_dict[tr_id] = gene_id
    else:
        tr_id = eachline.attr['Parent']
        gene_id = tr_gene_dict[tr_id]
    gene_name = gene_id
    new_attr = 'gene_id "{gene_id}"; transcript_id "{tr_id}"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "{gene_id}";'.format(**locals())
    gtf_file_info.write('%s\t%s\n' % (pos_info,new_attr))
gtf_file_info.close()



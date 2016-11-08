from Bio import SeqIO
import sys
import os
import re

if not len(sys.argv) == 3:
    print '    python ' + sys.argv[0] + ' pep.fasta annotation.inf'
    sys.exit(1)

pep_fasta = sys.argv[1]
annotation_out = sys.argv[2]

annotation_out_inf = open(annotation_out, 'w')
annotation_out_inf.write('Protein_ID\tTranscript_ID\tGene_ID\tGene_symbol\tGene_description\n')
for record in SeqIO.parse(pep_fasta, 'fasta'):
    description = record.description
    gene_id, transcript_id, gene_symbol, gene_des = ['--']*4
    if re.search(r'gene:(\S+)',description):
        gene_id = re.search(r'gene:(\S+)',description).groups()[0]
    if re.search(r'transcript:(\S+)',description):
        transcript_id = re.search(r'transcript:(\S+)',description).groups()[0]
    if re.search(r'gene_symbol:(\S+)',description):
        gene_symbol = re.search(r'gene_symbol:(\S+)',description).groups()[0]
    if re.search(r'description:(\S+)',description):
        gene_des = re.search(r'description:(.*)',description).groups()[0]
    annotation_out_inf.write('{record.id}\t{transcript_id}\t{gene_id}\t{gene_symbol}\t{gene_des}\n'.format(**locals()))
annotation_out_inf.close()
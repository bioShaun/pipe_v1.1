import sys
import os
from Bio import SeqIO

if not len(sys.argv) == 6:
    print '    python ' + sys.argv[0] + ' sp miRNA.seq utr.seq targetscan.input targetscan.utr'
    sys.exit(0)

sp = sys.argv[1]
miRNA_seq = sys.argv[2]
utr_seq = sys.argv[3]
targetscan_input = sys.argv[4]
targetscan_utr = sys.argv[5]

targetscan_input_inf = open(targetscan_input, 'w')
for record in SeqIO.parse(miRNA_seq, 'fasta'):
    miRNA_id = record.id
    miRNA_seq = str(record.seq)
    sp_inf = sp
    seed_region = miRNA_seq[1:8]
    targetscan_input_inf.write('%s\t%s\t%s\n' % (miRNA_id, seed_region, sp_inf))
targetscan_input_inf.close()

targetscan_utr_inf = open(targetscan_utr, 'w')
for record in SeqIO.parse(utr_seq, 'fasta'):
    utr_id = record.id
    utr_seq = str(record.seq)
    sp_inf = sp
    targetscan_utr_inf.write('%s\t%s\t%s\n' % (utr_id, sp_inf, utr_seq))
targetscan_utr_inf.close()

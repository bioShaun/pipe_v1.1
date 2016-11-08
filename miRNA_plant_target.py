'''
predict plant miRNA target by targetfinder

example:
targetfinder.pl -s uguuauuguggaaguucuuga -d /home/public/database/others/Cucumis_sativus/cucumber_ChineseLong_v2.cds.fa -p table -q test1

'''

import sys
import python_tools
from Bio import SeqIO
import os

if not len(sys.argv) == 4:
    print 'python ' + sys.argv[0] + ' miRNA_mature mRNA_sequence output_dir'
    sys.exit(0)

miRNA_mature = sys.argv[1]
mRNA_seq = sys.argv[2]
output_dir = sys.argv[3]

output_file = os.path.join(output_dir, 'miRNA_targetfinder_results.txt')
cmd_file = os.path.join(output_dir, 'miRNA_target_cmd.sh')

cmd_file_info  = open(cmd_file, 'w')
for seq_record in SeqIO.parse(miRNA_mature, "fasta"):
    seq_name = seq_record.id
    seq = str(seq_record.seq)
    cmd_file_info.write('targetfinder.pl -s {seq} -d {mRNA_seq} -p table -q {seq_name}\n'.format(**locals()))
cmd_file_info.close()

cmd = 'nohup sh {cmd_file} 1>{output_file} 2>{cmd_file}.log &'.format(**locals())
python_tools.circ_call_process(cmd)

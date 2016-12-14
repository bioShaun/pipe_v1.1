import sys
import os

if not len(sys.argv) == 3:
    print '    python ' + sys.argv[0] + ' fq.dir out.dir'
    sys.exit(0)

fq_dir = sys.argv[1]
out_dir = sys.argv[2]

clean_dir = os.path.join(fq_dir, 'Cleandata')
if not os.path.isdir(clean_dir):
    sys.exit('not find clean dir')
fq_dir_list = os.listdir(clean_dir)

for each_dir in fq_dir_list:
    each_sample_fq = os.path.join(clean_dir,'{0}/{0}_1.fq.gz'.format(each_dir))
    out_fq = os.path.join(out_dir, '{0}.fq'.format(each_dir))
    out_fa = os.path.join(out_dir, '{0}.fa'.format(each_dir))
    os.system('zcat {0} | head -40000 > {1}'.format(each_sample_fq, out_fq))
    os.system('fastq_to_fasta -i {0} -o {1}'.format(out_fq, out_fa))
    os.system('rm {0}'.format(out_fq))

import sys
import os
from HTSeq import GFF_Reader

stringtie_gtf = sys.argv[1]
output = sys.argv[2]

out_dir = os.path.split(output)[0]
if not os.path.exists(out_dir):
    os.system('mkdir -p %s ' % out_dir)

output_info = open(output,'w')
tr_fpkm_dict = {}
for eachline in GFF_Reader(stringtie_gtf):
    if eachline.type == 'transcript':
        tr_id = eachline.attr['transcript_id']
        tr_fpkm  = eachline.attr['FPKM']
        tr_fpkm_dict[tr_id] = tr_fpkm
        output = "%s;\n" % eachline.get_gff_line().strip()
    else:
        tr_id = eachline.attr['transcript_id']
        tr_fpkm = tr_fpkm_dict[tr_id]
        output = '%s; FPKM "%s";\n' % (eachline.get_gff_line().strip(), tr_fpkm)
    output_info.write(output)
output_info.close()

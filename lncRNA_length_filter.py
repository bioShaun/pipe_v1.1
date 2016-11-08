import sys
import os
import RNAseq_tools
import python_tools

if not len(sys.argv) == 3:
    print 'python ' + sys.argv[0] + ' gtf length_filtered_gtf'
    sys.exit(0)

ori_gtf = sys.argv[1]
out_gtf = sys.argv[2]

tr_dict = RNAseq_tools.get_transcript_info(ori_gtf)
filter_gtf_list = RNAseq_tools.gtf_length_filter(ori_gtf, tr_dict)
python_tools.write_obj_to_file(filter_gtf_list, out_gtf)
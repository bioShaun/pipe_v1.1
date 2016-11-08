import sys
import os
import re

if not len(sys.argv) == 3:
    print '    python ' + sys.argv[0] + ' geneid_output plot_data'
    sys.exit(0)

input_file = sys.argv[1]
output_file = sys.argv[2]

output_file_info = open(output_file,'w')
output_file_info.write('ID\tCoding_Potential\n')
with open(input_file) as input_file_info:
    for eachline in input_file_info:
        eachline_info = eachline.strip().split('\t')
        if len(eachline_info) == 8:
            each_id = eachline_info[-1].split()[-1]
            each_score = eachline_info[-3].strip()
            output_file_info.write('%s\t%s\n' % (each_id,each_score))
output_file_info.close()

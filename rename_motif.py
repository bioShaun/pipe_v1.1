import sys
import python_tools

if not len(sys.argv) == 3:
    print '    python ' + sys.argv[0] + ' fimo_out rename_fimo_out'
    sys.exit(0)

name_map = '/home/public/database/meme/pfm_plants.name'
name_map_dict = python_tools.table_to_dict(name_map,1,2,False)

fimo_out = sys.argv[1]
rename_fimo_out = sys.argv[2]
rename_fimo_out_info = open(rename_fimo_out, 'w')
with open(fimo_out) as fimo_out_info:
    for n, eachline in enumerate(fimo_out_info):
        if n == 0:
            rename_fimo_out_info.write(eachline)
        else:
            eachline_info = eachline.strip().split('\t')
            motif_id = eachline_info[0]
            motif_id_info = '\t'.join(eachline_info[1:])
            if motif_id in name_map_dict:
                motif_id = name_map_dict[motif_id]
            rename_fimo_out_info.write('%s\t%s\n' % (motif_id, motif_id_info)) 
rename_fimo_out_info.close()

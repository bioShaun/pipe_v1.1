'''
Usage:
extract_ncRNA_from_assemblyline.py <all.gtf> <ncRNA.gtf>

Extract ncRNA gtf from assemblyline assembled gtf

'''

from docopt import docopt
from HTSeq import GFF_Reader

if __name__ == '__main__':
    arguments = docopt(__doc__, version = "v1")
    all_gtf = arguments['<all.gtf>']
    nc_gtf = arguments['<ncRNA.gtf>']
    nc_gtf_inf = open(nc_gtf, 'w')
    for eachline in GFF_Reader(all_gtf):
        if 'transcript_category' in eachline.attr:
            if eachline.attr['transcript_category'] == 'ncRNA' or eachline.attr['transcript_category'] == 'lncRNA':
                nc_gtf_inf.write(eachline.get_gff_line())


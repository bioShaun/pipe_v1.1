import sys
import os
import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument('--kegg_table', help = 'KOBAS kegg enrichment analysis result', required = True)
parser.add_argument('--pathway_dir', help = 'pathway directory', required = True)
parser.add_argument('--log_file', help = 'check result', required = True)
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG,
                format='%(asctime)s %(filename)s [line:%(lineno)d] %(levelname)s : %(message)s',
                datefmt='%a, %d %b %Y %H:%M:%S',
                filename='%s' % args.log_file,
                filemode='w')


with open(args.kegg_table) as kegg_table_info :
    for n,eachline in enumerate(kegg_table_info) :
        if n != 0 :
            eachline_info = eachline.strip().split('\t')
            ko_id = eachline_info[2]
            pathway_pic1 = os.path.join(args.pathway_dir,'%s.pathview.png' % ko_id)
            pathway_pic2 = os.path.join(args.pathway_dir,'%s.pathview.gene.png' % ko_id)
            pathway_files = [pathway_pic1,pathway_pic2]
            for each in pathway_files :
                if not os.path.isfile(each) :
                    logging.error('%s not exists !' % each)
                elif not os.stat(each) :
                    logging.error('%s is empty !' % each)


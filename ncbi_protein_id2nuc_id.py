from Bio import Entrez
import sys
import os
import argparse
import json
import time
from urllib2 import HTTPError
import python_tools

cwd = os.getcwd()

parser = argparse.ArgumentParser(description="genome fasta stastics")
parser.add_argument('--protein',help='protein id file',required=True)
parser.add_argument('--protein_col',help='protein id column', default = '1')
parser.add_argument('--nucgi', help = 'nucleotide gi', required=True)
args = parser.parse_args()

Entrez.email = 'A.N.Other@example.com'

def accession2gi(accession) :
    try:
        handle = Entrez.elink(dbfrom="protein", id=accession, linkname="protein_nuccore")
    except HTTPError:
        time.sleep(20)
        handle = Entrez.elink(dbfrom="protein", id=accession, linkname="protein_nuccore")
    record = Entrez.read(handle)
    handle.close()
    linked = record[0]['LinkSetDb'][0]['Link'][0]['Id']
    time.sleep(10)
    return linked


if __name__ == '__main__':
    col_num = int(args.protein_col)
    gi_dict = {}
    with open(args.protein) as protein_info:
        for eachline in protein_info:
            eachline = eachline.strip()
            eachline_info = eachline.split('\t')
            pro_id = eachline_info[col_num-1]
            gi_dict[pro_id] = [eachline]

    for each_id in gi_dict:
        each_accession= "NA"
        try:
            each_accession = accession2gi(each_id)
        except Exception as e:
            print e
            print 'can not find gi number of accession: %s' % each_id
        gi_dict[each_id].append(each_accession)
    
    with open(args.nucgi,'w') as nucgi_info:
        for each_id in gi_dict:
            nucgi_info.write('%s\t%s\n' % (gi_dict[each_id][0], gi_dict[each_id][1]))



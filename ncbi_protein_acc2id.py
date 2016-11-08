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
parser.add_argument('--protein_acc',help='protein accession',required=True)
parser.add_argument('--protein_id', help = 'protein gi', required = True)
args = parser.parse_args()

Entrez.email = 'A.N.Other@example.com'

def accession2gi(accession) :
    try:
        handle = Entrez.esearch(db="protein", term = accession)
    except HTTPError:
        time.sleep(20)
        handle = Entrez.esearch(db="protein", term = accession)
    record = Entrez.read(handle)
    handle.close()
    linked = record['IdList'][0]
    time.sleep(10)
    return linked


if __name__ == '__main__':
    accession_list = [each.strip() for each in open(args.protein_acc)]
    gi_dict = {}
    for each_id in accession_list:
        each_accession= "NA"
        try:
            each_accession = accession2gi(each_id)
#            gi_list.append(accession2gi(each_id))
        except Exception as e:
            print e
            print 'can not find gi number of accession: %s' % each_id            
        gi_dict[each_id] = each_accession
    python_tools.write_obj_to_file(gi_dict,args.protein_id)

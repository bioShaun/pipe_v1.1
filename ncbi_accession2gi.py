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
parser.add_argument('--accession',help='nucleotide accession',required=True)
parser.add_argument('--gi', help = 'nucleotide gi', default = cwd)
args = parser.parse_args()

Entrez.email = 'A.N.Other@example.com'

def accession2gi(accession) :
    try:
        handle = Entrez.elink(dbfrom="nucleotide", id=accession, linkname="nucleotide_nuccore")
    except HTTPError:
        time.sleep(20)
        handle = Entrez.elink(dbfrom="nucleotide", id=accession, linkname="nucleotide_nuccore")
    record = Entrez.read(handle)
    handle.close()
    linked = record[0]['IdList'][0]
    time.sleep(10)
    return linked


if __name__ == '__main__':
    accession_list = [each.strip() for each in open(args.accession)]
    gi_list = []
    for each_id in accession_list:
        each_accession= ""
        try:
            each_accession = accession2gi(each_id)
#            gi_list.append(accession2gi(each_id))
        except:
            print 'can not find gi number of accession: %s' % each_id
        else:
            gi_list.append(each_accession)
    python_tools.write_obj_to_file(gi_list,args.gi)



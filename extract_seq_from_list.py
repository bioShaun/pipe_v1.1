import sys
import os
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--all_fa', help = 'all fasta sequence.',required =True)
parser.add_argument('--id_list', help = 'fasta id list.', required =True)
parser.add_argument('--ex', action='store_true' ,help = 'exclude fasta in id list.', default = False)
parser.add_argument('--sub_fa', help = 'exclude or extract sequence from all fasta file by id list.', required =True)
args = parser.parse_args()

id_list = [each.strip() for each in open(args.id_list,'r')]

my_records = []
for seq_record in SeqIO.parse(args.all_fa, "fasta") :
    seq_id = seq_record.id
    if args.ex :
        if seq_id not in id_list :
            my_records.append(seq_record)
    else :
        if seq_id in id_list :
            my_records.append(seq_record)
SeqIO.write(my_records, args.sub_fa, "fasta")


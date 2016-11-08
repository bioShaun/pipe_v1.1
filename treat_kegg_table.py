import RNAseq_tools
import sys
import os
import argparse

parser = argparse.ArgumentParser(description = 'Remove unnessecary information in KEGG table.')
parser.add_argument('--kegg_out', help = 'KOBAS KEGG analysis output table', required = True)
args = parser.parse_args()

my_kegg = RNAseq_tools.KEGG_enrich()
my_kegg.treat_KEGG_table(args.kegg_out)



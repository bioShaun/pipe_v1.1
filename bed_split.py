import argparse
import os

root_dir=os.getcwd()
##################################################################################
#parse the arguments
parser = argparse.ArgumentParser(description="genome fasta stastics")
parser.add_argument('--bed',help='bedfile',required=True)
parser.add_argument('--outdir',help='subbed directory',default=root_dir)

argv=vars(parser.parse_args())
bedfile=argv['bed'].strip()
outdir=argv['outdir'].strip()

bed_dict = {}

with open(bedfile,"r") as bed_info :
    for eachline in bed_info :
        eachline =eachline.strip()
        if eachline != "" :
            each_bed = eachline.split("\t")
            if each_bed[0] not in bed_dict :
                bed_dict[each_bed[0]] = []
                bed_dict[each_bed[0]].append(eachline)
            else :
                bed_dict[each_bed[0]].append(eachline)

for each_chr in bed_dict :
    outputfile = outdir+"/"+each_chr+".bed"
    with open(outputfile,"w") as chr_bed :
        for each_tr in bed_dict[each_chr] :
            chr_bed.write(each_tr+"\n")


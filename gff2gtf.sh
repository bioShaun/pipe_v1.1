gff=$1
gtf=$2

gffread $gff -F -O -E -T -o- > $gtf

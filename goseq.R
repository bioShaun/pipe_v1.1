library(goseq)
args<-commandArgs(T)

diff_gene <-  scan(args[1], what=character())
gene_length <- read.table(args[2],header=F,sep='\t')
go_id <- read.table(args[3],header=T,sep=',')
outfile<-args[4]

all_id<-gene_length[,1]
gene.vector = as.integer(all_id %in% diff_gene)
names(gene.vector) = all_id

id_len<-gene_length[,2]
names(id_len) = all_id

pwf=nullp(gene.vector,bias.data=id_len)
GO.wall=goseq(pwf,gene2cat=go_id)
GO.wall<-GO.wall[GO.wall$numDEInCat > 0,c(1,2,4,5,6,7)]
GO.wall$qvalue<-p.adjust(GO.wall$over_represented_pvalue,method="BH",n=length(GO.wall$over_represented_pvalue))
out_go<-GO.wall[,c(1,2,7,3,4,5,6)]
out_go<-na.omit(out_go)
if (dim(out_go)[1] > 0) {
    write.table(out_go,file=outfile,quote=F,sep='\t',row.names=F)
}else {
    print("No gene successfully annotated!")
}

# save(list=ls(all=TRUE), file='test.Rdata')

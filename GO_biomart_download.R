library("biomaRt")

args<-commandArgs(T)
gene_tr_file <- args[1]
data_name <- args[2]
out_name <- args[3]

geneid_info <- read.table(gene_tr_file, header = F, sep = '\t')
gene_id <- unique(geneid_info$V1)
ensembl = useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset(data_name, mart = ensembl)
goids = getBM(attributes=c('ensembl_gene_id','go_id'), filters='ensembl_gene_id', values=gene_id, mart=ensembl)
write.table(goids, file = paste(out_name, ".go.txt", sep = ""), quote = F, sep = ',' , row.names = F )

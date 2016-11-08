library("biomaRt")
ensembl = useMart("ensembl",dataset="ggallus_gene_ensembl")
gene_id<-"/media/zxchen/8679a838-8984-4175-82cd-4addacce31c5/project/OM-mRNA-chicken-P20160731/gene.list"
gene_id_inf<-read.table(gene_id, header = F, sep = '\t')
gene_inf<-getBM(attributes=c('ensembl_gene_id','external_gene_name','interpro'), filters='ensembl_gene_id', values=gene_id_inf, mart=ensembl)
write.table(gene_inf, file = paste("gga", ".des.tmp", sep = ""), quote = F, sep = ',' , row.names = F )
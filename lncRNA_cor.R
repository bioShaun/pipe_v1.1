#2016-8-4
library(dplyr)
library(WGCNA)
library(multtest)
#read data:
# data <- read.table("mRNA.exp.txt",header = T)
# #for test:
# data1 <- sample_n(data,100)
# data2 <- sample_n(data,100)
# data1 <- data1[,1:20]
# data2 <- data2[,1:20]
args<-commandArgs(T)
lncRNA_data <- args[1]
mRNA_data <- args[2]

data1 <- read.table(lncRNA_data,header = T, row.names = 1)
data2 <- read.table(mRNA_data,header = T, row.names = 1)
create_cor <- function(data1,data2){
  data1 <- as.matrix(data1)
  data2 <- as.matrix(data2)
  cors <- NULL;name1 <- NULL;name2 <- NULL
  for(i in 1:dim(data1)[1]){
    for(j in 1:dim(data2)[1]){
      cor <- cor(data1[i,],data2[j,])
      cors <- c(cors,cor)
      namei <- rownames(data1)[i]
      name1 <- c(name1,namei)
      namej <- rownames(data2)[j]
      name2 <- c(name2,namej)
    }
  }
  cor_data <- data.frame(lncRNA=name1,mRNA=name2,cor=cors)
}
cor_datas <- create_cor(data1,data2) %>% na.omit()

cor_datas$pvalue <- WGCNA::corPvalueFisher(cor_datas$cor,dim(data1)[2])
adj_pvalue_list <- multtest::mt.rawp2adjp(cor_datas$pvalue,proc = "Bonferroni")

adj_df <- as.data.frame(adj_pvalue_list$adjp)
adj_df$index <- adj_pvalue_list$index
cor_datas$adj_pvalue <- 1 
cor_datas$adj_pvalue[adj_df$index] <- adj_df$Bonferroni

#设置pcc pvalue阈值
pcc <- 0.05;pvalue <- 0.01
len <- dim(cor_datas)[1]
select_cor_datas <- cor_datas %>% arrange(desc(cor)) %>% 
  slice(c(1:floor(len*pcc),ceiling(len*(1-pcc)):len)) 
#%>%  dplyr::filter(adj_pvalue <= pvalue)
select_cor_datas <- select_cor_datas[which(select_cor_datas$adj_pvalue <= pvalue),]

write.table(cor_datas,'all.correlation.txt',sep = '\t', quote = F, row.names = F)
write.table(select_cor_datas,'significant.correlation.txt',sep = '\t', quote = F, row.names = F)



#2016-8-11 
library(dplyr)
library(WGCNA)
library(multtest)
Rcpp::sourceCpp('/home/lxgui/scripts/create_cor.cpp')

args<-commandArgs(T)
lncRNA_data <- args[1]
mRNA_data <- args[2]
#out_file <- args[3]

data1 <- read.table(lncRNA_data,header = T, row.names = 1)
data2 <- read.table(mRNA_data,header = T, row.names = 1)

#for test:
#data1 <- sample_n(data1,5100)
#data2 <- sample_n(data2,3000)
# data1 <- data1[,1:20]
# data2 <- data2[,1:20]
#优化代码
#create_cor_1 <- function(data1,data2){
#  data1 <- as.matrix(data1)
#  data2 <- as.matrix(data2)
#  #预先分配空间
#  cors <- rep(0,dim(data1)[1]*dim(data2)[1])
#  name1 <- name2 <- cors
#  for(i in 1:dim(data1)[1]){
#    for(j in 1:dim(data2)[1]){
#      cor.data <- cor(data1[i,],data2[j,])
#      cors[(i-1)*100+j] <- cor.data
#      name1[(i-1)*100+j] <- rownames(data1)[i]
#      name2[(i-1)*100+j] <- rownames(data2)[j]
#    }
#  }
#  cor_data <- data.frame(cor_name1=name1,cor_name2=name2,cor=cors)
#}

#create_cor_2 <- function(data1,data2){
#  data1 <- as.matrix(data1)
#  data2 <- as.matrix(data2)
#  #预先分配空间
#  cors <- rep(0,dim(data1)[1]*dim(data2)[1])
#  name1 <- name2 <- cors
#  if(dim(data1)[1] > dim(data2)[1]){
#    data.big <- data1
#    data.samll <- data2
#  }else{
#    data.big <- data2
#    data.samll <- data1
#  }
#  len <- dim(data.big)[1]
#  for(i in 1:dim(data.samll)[1]){
#    tmp <- apply(data.big,1,cor,y=data.samll[i,])
#    cors[(i*len-len+1):(i*len)] <- tmp
#    name1[(i*len-len+1):(i*len)] <- names(tmp)
#    name2[(i*len-len+1):(i*len)] <- rep(rownames(data.samll)[i],len)
#  }
#  cor_data <- data.frame(cor_name1=name1,cor_name2=name2,cor=cors)
#}

#create_cor_3 <- function(data1,data2){
#  data1 <- as.matrix(data1)
#  data2 <- as.matrix(data2)
#  #预先分配空间
#  cors <- rep(0,dim(data1)[1]*dim(data2)[1])
#  name1 <- name2 <- cors
#  len <- dim(data1)[1]
#  for(i in 1:dim(data1)[1]){
#    for(j in 1:dim(data2)[1]){
#      cor.data <- corC(data1[i,],data2[j,])
#      cors[(i-1)*len+j] <- cor.data
#      name1[(i-1)*len+j] <- rownames(data1)[i]
#      name2[(i-1)*len+j] <- rownames(data2)[j]
#    }
#  }
#  cor_data <- data.frame(cor_name1=name1,cor_name2=name2,cor=cors)
#}

#df1 <- create_cor_1(data1,data2)
data1.mat <- as.matrix(data1)
data2.mat <- as.matrix(data2)
cor_datas <- create_cor(data1.mat,data2.mat,rownames(data1),rownames(data2))
#cor_datas <- create_cor_3(data1,data2)
colnames(cor_datas) <- c('lncRNA', 'mRNA', 'cor')
#cor_datas <- na.omit(cor_datas)
#cor_datas <- cor_datas[cor_datas$lncRNA != cor_datas$mRNA,]
cor_datas$pvalue <- WGCNA::corPvalueFisher(round(cor_datas$cor,digits=5),dim(data1)[2])
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

write.table(cor_datas, file = 'all.correlation.txt', sep = '\t', quote = F, row.names = F, col.names = F)
write.table(select_cor_datas,'significant.correlation.txt',sep = '\t', quote = F, row.names = F)

#compare_df <- cbind(df1,df)
#code 瓶颈测试
#Rprof()
#create_cor(data1.mat,data2.mat,rownames(data1),rownames(data2)) %>% write.table("cor.txt",sep = "\t",quote = F,row.names = F)
# Rprof(NULL)
#summaryRprof()




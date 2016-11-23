#2016-10-18
#PCA plot
#####
library(ggthemes)
if(!require(FactoMineR)){
  install.packages('FactoMineR')
}
if(!require(factoextra)){
  install.packages('factoextra')
}
argv <- commandArgs(T)
PCA_plot_path <- argv[1]
PCA_group_path <- argv[2]
output_path <- argv[3]
#PCA_plot_path <- "/media/zxchen/3dddd6c4-2700-41a5-a677-15b165fa4e64/project/test_proj/mRNA_pipeline/version1/proj1/OM-mRNA-5-Chicken-Analysis/Analysis_results/quantification"
#list.files(PCA_plot_path)
PCA_file <- "genes.TMM.EXPR.matrix"
#PCA_group_file <- "group_sample"
PCA_data <- read.delim(PCA_plot_path,header = T,sep = "\t")
PCA_group <- read.delim(PCA_group_path,header = F,sep = "\t")
PCA_data <- t(PCA_data[,2:dim(PCA_data)[2]])
PCA_data_df <- as.data.frame(PCA_data)
group_data <- PCA_group$V1[match(PCA_group$V2,rownames(PCA_data_df))]
res.pca <- PCA(PCA_data_df,graph = F)
p <- fviz_pca_ind(res.pca,geom = 'point',addEllipses = T,
                  habillage = group_data,
                  ellipse.level = 0.9)+theme_calc()
ggsave(filename = paste(output_path,'PCA_plot.pdf',sep = '/'),plot = p,width = 8,height = 6)
ggsave(filename = paste(output_path,'PCA_plot.png',sep = '/'),plot = p,type="cairo-png",width = 8,height = 6)




#2016-10-18
#for mRNA report reads qulity barplot
#####
library(scales)
library(ggplot2)
library(ggthemes) 
argv <- commandArgs(T)
file_path <- argv[1]
output_path <- argv[2]
options(stringsAsFactors = F)
#file_path <- "/media/zxchen/3dddd6c4-2700-41a5-a677-15b165fa4e64/project/Analysis_project_backup/OM-mRNA-Medicago_truncatula-P20161206/OM-mRNA-15-Medicago_truncatula/mRNA_analysis_results/qc/reads_quality"
#output_path <- "./"
all_files <- list.files(file_path)
all_quality_data <- read.delim(paste(file_path,all_files[grep('^all.+txt$',all_files)],sep = "/"),header = T,sep = "\t")
qulity_data_files <- all_files[grep('*_reads_quality.txt',all_files)] 
#if(length(qulity_data_files) > 9){  
#  qulity_data_files <- qulity_data_files[1:9]
#}

qulity_data <- list()

BarPlot <- function(data,Title){
  col <- data$color
  names(col) <- col
  p<- ggplot(data,aes(x=Quality,y=Proportion,fill=color))+
    geom_bar(stat = 'identity')+scale_fill_manual(values = col)+
    geom_segment(aes(x=30,y=0,xend=30,yend=max(data$Proportion)),
                 colour='red',linetype='dashed',size=1)+theme_bw()+guides(fill = F)+
    scale_y_continuous(breaks = seq(from = 0,to = max(data$Proportion),by = 0.1),
  labels = scales::percent(seq(from = 0,to = max(data$Proportion),by = 0.1)))+
    xlab('Quality Score')+ggtitle(label = Title)
  p
}

for(i in 1:length(qulity_data_files)){
  qulity_data[[i]] <- read.delim(paste(file_path,qulity_data_files[i],sep = "/"))
}
for(i in 1:length(qulity_data)){
  qulity_data[[i]]$color <- ifelse(qulity_data[[i]]$Quality <= 30,'dodgerblue','navy')
}
#---
split_str <- function(strings,Split){
  for(i in 1:nchar(strings)){
    if(substr(strings,i,i) == Split){
      return(c(substr(strings,1,i-1),substr(strings,i+1,nchar(strings))))
    }
  }
}
for(i in 1:length(qulity_data)){
  plot_title <- split_str(qulity_data_files[i],Split = '.')[1]
  p <- BarPlot(qulity_data[[i]],Title = plot_title)
  ggsave(filename = paste(output_path,paste(plot_title,'barplot.png',sep = '_'),sep = "/"),type="cairo-png",plot = p,width = 8,height = 6)
  ggsave(filename = paste(output_path,paste(plot_title,'barplot.pdf',sep = '_'),sep = "/"),plot = p,width = 8,height = 6)
}
#----
if(length(unique(all_quality_data$Sample_ID)) >= 9){
  sample_id_set <- unique(all_quality_data$Sample_ID)[1:9]
} else {
  sample_id_set <- unique(all_quality_data$Sample_ID)
}
all_quality_data <- all_quality_data[(all_quality_data$Sample_ID %in% sample_id_set),]
all_quality_data$color <- ifelse(all_quality_data$Quality <= 30,'dodgerblue','navy')
col <- all_quality_data$color
names(col) <- col
all_quality_data_plot <- ggplot(all_quality_data,aes(x=Quality,y=Percent,fill=color))+
  geom_bar(stat = 'identity')+scale_fill_manual(values = col)+
  geom_segment(aes(x=30,y=0,xend=30,yend=max(all_quality_data$Percent)),
               colour='red',linetype='dashed',size=1)+theme_bw()+guides(fill = F)+
  scale_y_continuous(breaks = seq(from = 0,to = max(all_quality_data$Percent),by = 0.2),
                                                     labels = scales::percent(seq(from = 0,to = max(all_quality_data$Percent),by = 0.2)))+
  facet_wrap(~Sample_ID,ncol = 3)
ggsave(filename = paste(output_path,'all_quality_data_barplot.png',sep = '/'),type="cairo-png",plot = all_quality_data_plot,width = 10,height = 6)

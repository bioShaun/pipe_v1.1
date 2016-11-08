#2016-10-13
#plot1 boxplot
#plot2 violin plot
#plot3 density plot
#changed on 2016-11-07(print each one plot)
######
library(reshape2)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(gridExtra)
options(stringsAsFactors = F)
argv <- commandArgs(T)
Gene.tmp_file_path <- argv[1]
plot_output_path <- argv[2]
# Gene.tmp_file_path = 'Gene.tpm.xls'
# plot_output_path = './'


data <- read.delim(Gene.tmp_file_path,header = T)
# data <- data[,1:6]
data.m <- melt(data,id = 'id')
calc_color <- colorRampPalette(calc_pal()(12))(length(unique(data.m$variable)))
data.m$value <- log10(data.m$value + 1)

theme_set(theme_bw()+theme(axis.text.x = element_text(angle = -90,color = 'black',
                                                        vjust = 0.5,hjust = 0,
                                                        size = rel(length(calc_color)*0.02)),
                             axis.text.y = element_text(size = rel(length(calc_color)*0.02)),
                             legend.text = element_text(size = rel(length(calc_color)*0.02)),
                             legend.key = element_blank()))

boxplot1 <- ggplot(data.m,aes(x=variable,y=value,fill=variable))+
  geom_boxplot(notch = T)+guides(fill = F)+
  scale_fill_manual(values = calc_color)+xlab("")+ylab("")
boxplot2 <- ggplot(data.m,aes(x=variable,y=value,fill=variable))+
  geom_boxplot(notch = T)+guides(fill = guide_legend(nrow = 8,title = 'sample'))+
  scale_fill_manual(values = calc_color)+xlab("")+ylab("")
violin <- ggplot(data.m,aes(x=variable,y=value,fill=variable))+
  geom_violin()+guides(fill = guide_legend(nrow = 8,title = 'sample'))+
  scale_fill_manual(values = calc_color)+xlab("")+ylab("")
density_plot <- ggplot(data.m,aes(value,fill = variable))+
  geom_density(alpha = .4)+
  scale_fill_manual(values = calc_color)+
  theme(axis.text.x = element_text(angle = 0))+
  guides(fill = guide_legend(nrow = 8,title = 'sample'))+xlab("")
#----print out----
#----each plot----

width_parameter = 0.6
heigth_parameter = 1.05

ggsave(filename = paste(plot_output_path,'Gene_expression.boxplot.png',sep = '/'),type="cairo-png",plot = boxplot2,
       width = length(unique(data.m$variable)) * width_parameter,height = 10)
ggsave(filename = paste(plot_output_path,'Gene_expression.boxplot.pdf',sep = '/'),plot = boxplot2,
       width = length(unique(data.m$variable)) * width_parameter,height = 10)
ggsave(filename = paste(plot_output_path,'Gene_expression.violin.png',sep = '/'),type="cairo-png",plot = violin,
       width = length(unique(data.m$variable)) * width_parameter,height = 10)
ggsave(filename = paste(plot_output_path,'Gene_expression.violin.pdf',sep = '/'),plot = violin,
       width = length(unique(data.m$variable)) * width_parameter,height = 10)
ggsave(filename = paste(plot_output_path,'Gene_expression.density_plot.png',sep = '/'),type="cairo-png",plot = density_plot,
       width = length(unique(data.m$variable)) * width_parameter,height = 10)
ggsave(filename = paste(plot_output_path,'Gene_expression.density_plot.pdf',sep = '/'),plot = density_plot,
       width = length(unique(data.m$variable)) * width_parameter,height = 10)
#----merge plot----
p <- grid.arrange(boxplot1,violin,density_plot,nrow = 2,
             layout_matrix = rbind(c(1,2),c(3,3)))
ggsave(filename = 'Gene_expression.png',type="cairo-png",plot = p,
       width = length(unique(data.m$variable)) * width_parameter,height = 10)
ggsave(filename = 'Gene_expression.pdf',plot = p,
       width = length(unique(data.m$variable)) * width_parameter,height = 10)


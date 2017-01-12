#setwd("C:/Users/Administrator/Desktop/R/data/数据")
#cwd = getwd()
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
library(argparser, quietly = TRUE)

#----set params-----

options(bitmapType='cairo')
options(stringsAsFactors = F)

parser <- arg_parser('mRNA report Volcano plot')
parser <- add_argument(parser, '--diff', help = 'merged diff file')
parser <- add_argument(parser, '--out_dir', help = 'output directory')
parser <- add_argument(parser, '--qvalue', help = 'significant qvlaue threshold', default = 0.05)
parser <- add_argument(parser, '--logfc', help = 'significant logfc threshold', default = 1)
argv <- parse_args(parser)

#
test_data_path <- argv$diff
p = argv$qvalue
lgfc = argv$logfc
out_dir <- argv$out_dir

y_line_pos <- round(-log10(p),1)
merge_data <- fread(test_data_path)
dpa_results <- merge_data[,c("logFC","FDR","Compare"),with =F]
dpa_results$logFDR <- -log10(dpa_results$FDR)
dpa_results$color <- "blue"
dpa_results$color <- ifelse(dpa_results$logFC > lgfc & dpa_results$FDR < p,"red",
                            ifelse(dpa_results$logFC < -(lgfc) & dpa_results$FDR < p,"green","blue"))
#theme set
theme_set(theme_bw()+theme(
   axis.text = element_text(colour = "black"),
   axis.title = element_text(colour = "black"),
   legend.position = c(0.2,0.9),
   legend.background = element_blank(),
   legend.text = element_text(face = "bold",colour = "black",size = rel(1.5)),
   legend.title = element_text(face = "bold",colour = "black",size = rel(1.5))
 ))


# # labels and titles
# Volcano_plot <- function(dpa_result,name){
#   count <- as.data.frame(table(dpa_result$color))
#   red_count <- sum(count$Freq[which(count$Var1 == "red")])
#   green_count <- sum(count$Freq[which(count$Var1 == "green")])
#   all_count_number <- red_count+green_count
#   all_count <- paste("Different Expressed Genes",all_count_number,sep = ":")
#   if(all_count_number < 100) {
#     dpa_result$logFC  <- ifelse(dpa_result$logFC > 15,
#                                   15+(dpa_result$logFC - 15)/logFC_max,
#                                   dpa_result$logFC)
#     dpa_result$logFC  <- ifelse(dpa_result$logFC < -15,
#                                   -15-(dpa_result$logFC + 15)/logFC_min,
#                                   dpa_result$logFC)
#     logFDR_max = max(dpa_result$logFDR)
#     dpa_result$logFDR <- ifelse(dpa_result$logFDR > 50,
#                                   50+(dpa_result$logFDR - 50)/logFDR_max,
#                                   dpa_result$logFDR)
#   }
#   logFC_max = max(dpa_result$logFC)
#   logFC_min = min(dpa_result$logFC)
#   logFC_limit <- ceiling(max(c(abs(logFC_min),logFC_max)))
#   logFC_limit <- ifelse(logFC_limit<8,8,logFC_limit)
#   logFC_limit <- ifelse(logFC_limit>15,15,logFC_limit)
#   logFDR_limit <- ceiling(max(dpa_result$logFDR))
#   logFDR_limit <- ifelse(logFDR_limit>50,50,logFDR_limit)
#   logFDR_limit <- ifelse(logFDR_limit<35,35,logFDR_limit)
#   red_label <- paste("up regulated",red_count,sep = ":")
#   green_label <- paste("down regulated",green_count,sep = ":")
#   #brewer for color
#   red <- brewer.pal(6,"Reds")[6]
#   green <- brewer.pal(6,"Greens")[6]
#   blue <- brewer.pal(4,"Blues")[4]
#   p <- ggplot(dpa_result,aes(logFC,logFDR,colour = color))+geom_point(size = .6)+
#     scale_color_manual(values = c("red"=red,"green"=green,"blue"=blue),
#                        breaks = c("red","green"),
#                        labels = c(red_label,green_label))+
#     geom_hline(yintercept = y_line_pos,lty = 4,size = .45)+
#     geom_vline(xintercept = -(lgfc),lty = 4,size = .45)+
#     geom_vline(xintercept = lgfc,lty = 4,size = .45)+
#     xlab("logFC")+ylab("-log10(FDR)")+ggtitle(name)+
#     scale_y_continuous(breaks = c(0,1.3,3,10,20,30),limits = c(0,logFDR_limit))+
#     scale_x_continuous(breaks = c(-8,-4,-2,-1,0,1,2,4,8),limits = c(-logFC_limit,logFC_limit))+
#     labs(color = all_count)
#   p
# }
# compare_sets <- unique(dpa_results$Compare)
# p1 <- Volcano_plot(dpa_results[Compare == compare_sets[1]],name = compare_sets[1])
# p2 <- Volcano_plot(dpa_results[Compare == compare_sets[1]],name = compare_sets[1])
# p3 <- Volcano_plot(dpa_results[Compare == compare_sets[1]],name = compare_sets[1])
# p4 <- Volcano_plot(dpa_results[Compare == compare_sets[1]],name = compare_sets[1])
# p5 <- Volcano_plot(dpa_results[Compare == compare_sets[1]],name = compare_sets[1])
# p6 <- Volcano_plot(dpa_results[Compare == compare_sets[1]],name = compare_sets[1])
# p7 <- Volcano_plot(dpa_results[Compare == compare_sets[1]],name = compare_sets[1])
# p8 <- Volcano_plot(dpa_results[Compare == compare_sets[1]],name = compare_sets[1])
# p9 <- Volcano_plot(dpa_results[Compare == compare_sets[1]],name = compare_sets[1])
# merge_plot <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,labels=c(compare_sets[1],compare_sets[2],compare_sets[3],
#                                                                      compare_sets[4],compare_sets[5],compare_sets[6],
#                                                                      compare_sets[7],compare_sets[8],compare_sets[9]),ncol = 3, nrow = 3)
# 
# ggsave(filename=paste(out_dir,"merge_Volcano_plot.pdf",sep="/"),plot=merge_plot,width = 15,height = 10)
# ggsave(filename=paste(out_dir,"merge_Volcano_plot.png",sep="/"),type="cairo-png",plot=merge_plot,width = 15,height = 10)
#---plan B----

#dpa_result <- dpa_results[dpa_results$Compare %in% select_sets]
dpa_result <- dpa_results
count <- as.data.frame(table(dpa_result$color))
red_count <- sum(count$Freq[which(count$Var1 == "red")])
green_count <- sum(count$Freq[which(count$Var1 == "green")])
all_count_number <- red_count+green_count
all_count <- paste("Different Expressed Genes",all_count_number,sep = ":")
if(all_count_number < 100) {
  dpa_result$logFC  <- ifelse(dpa_result$logFC > 15,
                              15+(dpa_result$logFC - 15)/logFC_max,
                              dpa_result$logFC)
  dpa_result$logFC  <- ifelse(dpa_result$logFC < -15,
                              -15-(dpa_result$logFC + 15)/logFC_min,
                              dpa_result$logFC)
  logFDR_max = max(dpa_result$logFDR)
  dpa_result$logFDR <- ifelse(dpa_result$logFDR > 50,
                              50+(dpa_result$logFDR - 50)/logFDR_max,
                              dpa_result$logFDR)
}
logFC_max = max(dpa_result$logFC)
logFC_min = min(dpa_result$logFC)
logFC_limit <- ceiling(max(c(abs(logFC_min),logFC_max)))
logFC_limit <- ifelse(logFC_limit<8,8,logFC_limit)
logFC_limit <- ifelse(logFC_limit>15,15,logFC_limit)
logFDR_limit <- ceiling(max(dpa_result$logFDR))
logFDR_limit <- ifelse(logFDR_limit>50,50,logFDR_limit)
logFDR_limit <- ifelse(logFDR_limit<35,35,logFDR_limit)
red_label <- paste("up regulated",red_count,sep = ":")
green_label <- paste("down regulated",green_count,sep = ":")
#brewer for color
red <- brewer.pal(6,"Reds")[6]
green <- brewer.pal(6,"Greens")[6]
blue <- brewer.pal(4,"Blues")[4]
p <- ggplot(dpa_result,aes(logFC,logFDR,colour = color))+geom_point(size = .6)+
  scale_color_manual(values = c("red"=red,"green"=green,"blue"=blue),
                     breaks = c("red","green"),
                     labels = c(red_label,green_label))+
  geom_hline(yintercept = y_line_pos,lty = 4,size = .45)+
  geom_vline(xintercept = -(lgfc),lty = 4,size = .45)+
  geom_vline(xintercept = lgfc,lty = 4,size = .45)+
  xlab("logFC")+ylab("-log10(FDR)")+ggtitle("")+
  scale_y_continuous(breaks = c(1.3,10,20,30),limits = c(0,logFDR_limit))+
  scale_x_continuous(breaks = c(-8,-4,-1,0,1,4,8),limits = c(-logFC_limit,logFC_limit))+
  labs(color = all_count)+facet_wrap(~Compare)+guides(color=F)
ggsave(filename=paste(out_dir,"merge_Volcano_plot.pdf",sep="/"),plot=p,width = 8,height = 8)
ggsave(filename=paste(out_dir,"merge_Volcano_plot.png",sep="/"),type="cairo-png",plot=p,width = 8,height = 8,dpi = 300)









  

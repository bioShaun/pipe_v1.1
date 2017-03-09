#setwd("C:/Users/Administrator/Desktop/R/data/数据")
cwd = getwd()
library(ggplot2)
library(data.table)
library(RColorBrewer)
options(bitmapType='cairo')

args<-commandArgs(T)
diff_file <- args[1]
name <- args[2]
out_dir <- args[3]
p <- as.numeric(args[4])
lgfc <- as.numeric(args[5])
# diff_file <- 'C-b_vs_F-a.edgeR.DE_results.txt'
# name <- 'C-b_vs_F-a'
# out_dir <- cwd

y_line_pos = round(-log10(p),1)

diff_file_info <- fread(diff_file)
dpa_results <- diff_file_info[,c("logFC","FDR"),with =F]
dpa_results$logFDR <- -log10(dpa_results$FDR)
dpa_results$color <- "blue"

up_name = unlist(strsplit(name,split = "_vs_"))[2]
down_name = unlist(strsplit(name,split = "_vs_"))[1]

for(i in 1:dim(dpa_results)[1]){
  if(dpa_results$logFC[i] > lgfc & dpa_results$FDR[i] < p)
      dpa_results$color[i] <- "red"
  else if(dpa_results$logFC[i] < -(lgfc) & dpa_results$FDR[i] < p)
      dpa_results$color[i] <- "green"
}

# labels and titles
dpa_results2 <- dpa_results
count <- table(dpa_results2$color)
count <- as.data.frame(count)

red_count <- sum(count$Freq[which(count$Var1 == "red")])
green_count <- sum(count$Freq[which(count$Var1 == "green")])
all_count_number <- red_count+green_count
all_count <- paste("Different Expressed Genes",all_count_number,sep = ":")


if (all_count_number < 100) {
dpa_results2$logFC  <- ifelse(dpa_results2$logFC > 15,
                               15+(dpa_results2$logFC - 15)/logFC_max,
                              dpa_results2$logFC)
dpa_results2$logFC  <- ifelse(dpa_results2$logFC < -15,
                              -15-(dpa_results2$logFC + 15)/logFC_min,
                              dpa_results2$logFC)
logFDR_max = max(dpa_results2$logFDR)
dpa_results2$logFDR <- ifelse(dpa_results2$logFDR > 50,
                              50+(dpa_results2$logFDR - 50)/logFDR_max,
                              dpa_results2$logFDR)
}

logFC_max = max(dpa_results2$logFC)
logFC_min = min(dpa_results2$logFC)
logFC_limit <- ceiling(max(c(abs(logFC_min),logFC_max)))
logFC_limit <- ifelse(logFC_limit<8,8,logFC_limit)
logFC_limit <- ifelse(logFC_limit>15,15,logFC_limit)
logFDR_limit <- ceiling(max(dpa_results2$logFDR))
logFDR_limit <- ifelse(logFDR_limit>50,50,logFDR_limit)
logFDR_limit <- ifelse(logFDR_limit<35,35,logFDR_limit)

#red2_count <- count$Freq[which(count$Var1 == "red2")]
#green2_count <- count$Freq[which(count$Var1 == "green2")]
#red1_label <- paste("up regulated",red1_count,sep = ":")
#green1_label <- paste("down regulated",green1_count,sep = ":")
red_label <- paste("No.", up_name, "up-regulated genes:",red_count,sep = " ")
green_label <- paste("No.", down_name, "down regulated genes:",green_count,sep = " ")
#brewer for color
red <- brewer.pal(6,"Reds")[6]
green <- brewer.pal(6,"Greens")[6]
blue <- brewer.pal(4,"Blues")[4]
#red1 <- reds[4];red2 <- reds[6]
#green1 <- greens[4];green2 <- greens[6]
#theme set
theme_set(theme_bw())
#plot
p <- ggplot(dpa_results2,aes(logFC,logFDR,colour = color))+geom_point(size = .6)+
  scale_color_manual(values = c("red"=red,"green"=green,"blue"=blue),
                     breaks = c("red","green"),
                     labels = c(red_label,green_label))+
  geom_hline(yintercept = y_line_pos,lty = 4,size = .45)+
  geom_vline(xintercept = -(lgfc),lty = 4,size = .45)+
  geom_vline(xintercept = lgfc,lty = 4,size = .45)+
  xlab("logFC")+ylab("-log10(FDR)")+ggtitle(name)+
  scale_y_continuous(breaks = c(0,1.3,3,10,20,30),limits = c(0,logFDR_limit))+
  scale_x_continuous(breaks = c(-8,-4,-2,-1,0,1,2,4,8),limits = c(-logFC_limit,logFC_limit))+
  labs(color = all_count)



ggsave(filename=paste(out_dir,"/",name,".Volcano_plot.pdf",sep=""),plot=p,width = 8,height = 6)
ggsave(filename=paste(out_dir,"/",name,".Volcano_plot.png",sep=""),type="cairo-png",plot=p,width = 8,height = 6)

# ggsave("15dpa_vs_m-20dpa.pdf",plot = p,width = 8,height = 6)
# ggsave()
# ggplot(dpa_results,aes(logFC,FDR,colour = color))+
#   geom_point(size = .6)+scale_y_continuous(breaks = c(0,1.3,3,10,20),limits = c(0,40))+
#   scale_x_continuous(limits = c(-4,4),breaks = c(-4,-2,-1,0,1,2,4))











  

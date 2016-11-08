library(ggplot2)
library(ggthemes)
library(scales)

args<-commandArgs(T)
length_file <- args[1]
name <- args[2]

length_file <- read.delim(length_file, sep = "\t", header = T)
length_file$prop <- length_file$Count/sum(length_file$Count)

p <- ggplot(data=length_file, aes(x=X.Length, y=prop,fill = prop)) +
  geom_bar(stat="identity",width = 0.75,colour = "black")
p <- p + scale_x_continuous(limits=c(18,35),breaks=seq(18,35))
p <- p + scale_y_continuous(labels = percent)
p <- p + theme_calc()+scale_fill_gradient2(low='red',mid="white",high="dodgerblue")
p <- p + guides(fill = F) + xlab("Length") + ylab("Proportion") + ggtitle(paste(name,"Length Distribution", sep = " "))

ggsave(filename=paste(name, '_length_distribution.png', sep = ""), type = 'cairo-png', plot = p, width = 8, height = 6, dpi=300)
ggsave(filename=paste(name, '_length_distribution.pdf', sep = ""), plot = p, width = 8, height = 6)
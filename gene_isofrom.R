library(ggplot2)
library(RColorBrewer)
library(scales)
args<-commandArgs(T)
gene_isofrom <- args[1]
gene_stat <- read.table(gene_isofrom,header = T,sep='\t',row.names = 1)
#gene_stat <- read.table('gene_isofrom.txt',header = T,sep='\t',row.names = 1)

gene_stat_plot <- subset(gene_stat,Type == "protein_coding" | Type == "lncRNA")

tr_lnc_plot <- subset(gene_stat_plot,Type == "lncRNA")
tr_mrna_plot <- subset(gene_stat_plot,Type == "protein_coding")
df_lnc_table <- table(tr_lnc_plot)
df_mrna_table <- table(tr_mrna_plot)
mytable_lnc <- prop.table(df_lnc_table)
mytable_mrna <- prop.table(df_mrna_table)

df_lnc_mytable <- as.data.frame(mytable_lnc,stringsAsFactors = F)
df_mrna_mytable <- as.data.frame(mytable_mrna,stringsAsFactors = F)
df_mytable2 <- rbind(df_lnc_mytable,df_mrna_mytable)
df_mytable2 <- subset(df_mytable2,Type == "protein_coding" | Type == "lncRNA")
df_mytable3 <- df_mytable2
df_mytable3$Isoform_number <- as.numeric(df_mytable3$Isoform_number)

p <- ggplot(df_mytable3, aes(x = Isoform_number,y=Freq,fill = Type)) +  
  geom_bar(position="dodge",stat = "identity")
p <- p + scale_y_continuous(labels = percent)
p <- p + scale_x_continuous(limits=c(0,9),breaks = c(1:8)) 
p <- p + scale_fill_manual(values=c("#313695", "#A50026")) + guides(fill=guide_legend(title=NULL))
p <- p + xlab('Number of isoform(s)') + ylab('Proportion of genes(%)') + labs(title="Number of alternatively spliced isoform(s) per gene locus")
p
ggsave(filename='Number_of_alternatively_spliced_isoform.png',type = 'cairo-png',plot = p ,width = 8, height = 8,dpi=300)
ggsave(filename='Number_of_alternatively_spliced_isoform.pdf',plot = p ,width = 8, height = 8)

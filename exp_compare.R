library(ggplot2)
library(RColorBrewer)
library(scales)

args<-commandArgs(T)
mRNA_exp <- args[1]
lncRNA_exp <- args[2]

protein_coding <- read.table(mRNA_exp,header= T,row.names = 1, sep = "\t")
protein_coding_plot <- rowMeans(protein_coding)
protein_coding_plot <- as.data.frame(protein_coding_plot)
lg_protein_coding_plot <- protein_coding_plot
lg_protein_coding_plot$protein_coding_plot <- log10(lg_protein_coding_plot$protein_coding_plot)
lg_protein_coding_plot$protein_coding_plot <- na.omit(lg_protein_coding_plot$protein_coding_plot)
lg_protein_coding_plot$protein_coding_plot <- ceiling(lg_protein_coding_plot$protein_coding_plot)
lg_protein_coding_plot$bio_type <- "protein_coding"
colnames(lg_protein_coding_plot)<- c("value","bio_type")
df_mrna_table <- table(lg_protein_coding_plot)
mytable_mrna <- prop.table(df_mrna_table)
df_mrna_mytable <- as.data.frame(mytable_mrna,stringsAsFactors = F)

lnc <- read.table(lncRNA_exp,header= T,row.names = 1, sep = "\t")
lnc_plot <- rowMeans(lnc)
lnc_plot <- as.data.frame(lnc_plot)
lg_lnc_plot <- lnc_plot
lg_lnc_plot$lnc_plot <- log10(lg_lnc_plot$lnc_plot)
lg_lnc_plot$lnc_plot <- na.omit(lg_lnc_plot$lnc_plot )
lg_lnc_plot$lnc_plot <- ceiling(lg_lnc_plot$lnc_plot )
lg_lnc_plot$bio_type <- "lncRNA"
colnames(lg_lnc_plot)<- c("value","bio_type")

df_lnc_table <- table(lg_lnc_plot)
mytable_lnc <- prop.table(df_lnc_table)
df_lnc_mytable <- as.data.frame(mytable_lnc,stringsAsFactors = F)

df_plot <- rbind(df_mrna_mytable,df_lnc_mytable)

df_plot$value <- as.numeric(df_plot$value)
df_plot <- df_plot[which(df_plot$value != -Inf),]

p <- ggplot(df_plot, aes(x=value,y=Freq,fill = bio_type)) + geom_bar(stat = 'identity',position=position_dodge(0.85),width = 0.85)
p <- p + scale_x_continuous(limits=c(-4,4),breaks = c(-3:3)) 
p <- p + scale_y_continuous(labels = percent)
p <- p + scale_fill_manual(values=c("#313695", "#A50026")) + guides(fill=guide_legend(title=NULL))
p <- p + xlab('Expression bins (with 1 log10 units)') + ylab('Proportion of genes(%)') + labs(title="Proportion of lncRNA/mRNA by log10 TPM")

ggsave(filename='Proportion_of_lncRNA_mRNA_by_log10_TPM.png',type = 'cairo-png',plot = p ,width = 8, height = 8,dpi=300)
ggsave(filename='Proportion_of_lncRNA_mRNA_by_log10_TPM.pdf',plot = p ,width = 8, height = 8)
library(RColorBrewer)
library(ggplot2)
library(reshape)
args<-commandArgs(T)
gene_exp <- args[1]
out_dir <- args[2]

df <- read.table(file=gene_exp, row.names = 1, sep = '\t', header = T)
df_mat <- as.matrix.data.frame(df)
log_df_mat <- log10(df_mat+1)
df_plot <- as.data.frame(log_df_mat)
df_plot2 <- melt(df_plot)

p <- ggplot(df_plot2, aes(x=variable, y=value, fill = variable)) + 
  geom_boxplot(notch=TRUE) + theme_bw()
p <- p +scale_fill_brewer(palette="Dark2")+ guides(fill=guide_legend(title=NULL))
p <- p + xlab('') + ylab('log10(TPM+1)') + labs(title="") + theme_bw()
ggsave(filename=paste(out_dir, '/Gene_expression.png', sep = ''),type = 'cairo-png',plot = p ,width = 8, height = 6,dpi=300)
ggsave(filename=paste(out_dir, '/Gene_expression.pdf', sep = ''),plot = p ,width = 8, height = 6)
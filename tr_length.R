library(ggplot2)
library(RColorBrewer)
library(scales)
args<-commandArgs(T)

transcript_feature <- args[1]
tr_feature <- read.table(transcript_feature, header = T,sep='\t')

#tr_feature <- read.table('transcript_feature.txt',header = T,sep='\t')
tr_feature_plot <- subset(tr_feature,Type == "protein_coding" | Type == "lncRNA")
len_max <- round(quantile(tr_feature_plot$Length,0.98),-1)

## length distribution
p <- ggplot(tr_feature_plot,aes(Length,fill = Type)) + geom_density(alpha = 0.8) + xlim(0,len_max)
p <- p + scale_fill_manual(values=c("#313695", "#A50026"))+guides(fill=guide_legend(title=NULL))
p <- p + xlab('Size (bp)') + ylab('Density') + labs(title="Processed Transcript Size Distribution")
# p <- p + facet_wrap(~status)
# p
ggsave(filename='Processed_Transcript_Size_Distribution.png',type = 'cairo-png',plot = p ,width = 8, height = 8,dpi=300)
ggsave(filename='Processed_Transcript_Size_Distribution.pdf',plot = p ,width = 8, height = 8)

## exon number 
tr_exon_plot <- tr_feature_plot[,c(3,4)]
tr_lnc_plot <- subset(tr_exon_plot,Type == "lncRNA")
tr_mrna_plot <- subset(tr_exon_plot,Type == "protein_coding")
df_lnc_table <- table(tr_lnc_plot)
df_mrna_table <- table(tr_mrna_plot)
mytable_lnc <- prop.table(df_lnc_table)
mytable_mrna <- prop.table(df_mrna_table)


df_lnc_mytable <- as.data.frame(mytable_lnc,stringsAsFactors = F)
df_mrna_mytable <- as.data.frame(mytable_mrna,stringsAsFactors = F)
df_mytable2 <- rbind(df_lnc_mytable,df_mrna_mytable)
df_mytable2 <- subset(df_mytable2,Type == "protein_coding" | Type == "lncRNA")
df_mytable3 <- df_mytable2
df_mytable3$Exon_number <- as.numeric(df_mytable3$Exon_number)


p <- ggplot(df_mytable3, aes(x = Exon_number,y=Freq,fill = Type)) +  
     geom_bar(position="dodge",stat = "identity")
p <- p + scale_y_continuous(labels = percent)
p <- p + scale_x_continuous(limits=c(0,16),breaks = c(1:15)) 
p <- p + scale_fill_manual(values=c("#313695", "#A50026")) + guides(fill=guide_legend(title=NULL))
p <- p + xlab('Number of exon(s)') + ylab('Proportion of transcripts(%)') + labs(title="Exon per transcript")
p
ggsave(filename='Exon_per_transcript.png',type = 'cairo-png',plot = p ,width = 8, height = 8,dpi=300)
ggsave(filename='Exon_per_transcript.pdf',plot = p ,width = 8, height = 8)

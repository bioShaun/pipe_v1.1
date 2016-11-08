library(ggplot2)
library(RColorBrewer)
args<-commandArgs(T)
exon_data <- args[1]
intron_data <- args[2]
lnc_exon <- read.table(exon_data, header = T,sep='\t')
lnc_intron <- read.table(intron_data, header = T,sep='\t')

# lnc_exon <- read.table('exon_length.txt',header = T,sep='\t')
# lnc_intron <- read.table('intron_length.txt',header = T,sep='\t')
lnc_exon$status <- "exon"
lnc_intron$status <- "intron"
lnc_exon_intron <- rbind(lnc_exon,lnc_intron)
lnc_exon_intron_plot <- subset(lnc_exon_intron,Type == "protein_coding" | Type == "lncRNA")
len_max <- round(quantile(lnc_exon_intron_plot$Length,0.98),-1)

mypal <- colorRampPalette( brewer.pal( 3 , "RdBu" ) )

p <- ggplot(lnc_exon_intron_plot,aes(Length,fill = Type)) + geom_density(alpha = 0.8) + xlim(0,2000)
p <- p + scale_fill_manual(values=c("#313695", "#A50026"))+guides(fill=guide_legend(title=NULL))
p <- p + xlab('Size(bp)') + ylab('Density') + labs(title="Exon/Intron Size Distribution")
p <- p + facet_wrap(~status)
#p
ggsave(filename='Exon_Intron_Size_Distribution.png',type = 'cairo-png',plot = p ,width = 8, height = 8,dpi=300)
ggsave(filename='Exon_Intron_Size_Distribution.pdf',plot = p ,width = 8, height = 8)

library(ggplot2)
library(RColorBrewer)
library(scales)
#library("ddply")


args<-commandArgs(T)
mRNA_cp_file = args[1]
lnc_cp_file = args[2]

mRNA_cp <- read.table(mRNA_cp_file,header= T, sep = "\t")
lnc_cp <- read.table(lnc_cp_file,header= T, sep = "\t")

# mRNA_cp <- read.table("mRNA_candidat.codingpred.plot",header= T, sep = "\t")
# lnc_cp <- read.table("lncRNA_candidat.codingpred.plot",header= T, sep = "\t")
lnc_cp$bio_type <- "lncRNA"
mRNA_cp$bio_type <- "protein_coding"
merged_df <- rbind(mRNA_cp, lnc_cp)

p <- ggplot(merged_df, aes(x=bio_type,y=Coding_Potential,fill = bio_type))
p <- p+ geom_boxplot(outlier.shape=NA, notch=TRUE) + ylim(-10,160)
p <- p + scale_fill_manual(values=c("#313695", "#A50026")) + guides(fill=guide_legend(title=NULL))
p <- p + xlab('') + ylab('GeneId Coding Potential') + labs(title="") + theme_bw()

ggsave(filename='GeneId_Coding_Potential.png',type = 'cairo-png',plot = p ,width = 6, height = 8,dpi=300)
ggsave(filename='GeneId_Coding_Potential.pdf',plot = p ,width = 6, height = 8)

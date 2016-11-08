args<-commandArgs(TRUE)
name = args[1]
GOdata.filepath = args[2]
KEGGdata.filepath = args[3]
outdir = args[4]
library(ggplot2)
#load data
#GOdata.filepath <- "C:/Users/Administrator/Desktop/R/data/数据/6023-0dpi_vs_6023-1dpi.GO.enrich.xls"
#KEGGdata.filepath <- "C:/Users/Administrator/Desktop/R/data/数据/6023-0dpi_vs_6023-1dpi.KEGG.enrich.xls"
GO <- read.delim(GOdata.filepath,head=T,stringsAsFactors=F)
KEGG <- read.delim(KEGGdata.filepath,head=T,stringsAsFactors=F)
names(KEGG)[1] <- "Name"
names(KEGG)[7] <- "P.Value"
names(GO)[3] <- "P.Value"
names(GO)[6] <- "Name"

for(i in 1:dim(KEGG)[1]){
  if(nchar(KEGG$Name[i]) > 50){
    KEGG$Name[i] <- substr(KEGG$Name[i],1,50)
  }
}
for(i in 1:dim(GO)[1]){
  if(nchar(GO$Name[i]) > 50){
    GO$Name[i] <- substr(GO$Name[i],1,50)
  }
}

#----set figure style----
# old<-theme_bw()+theme(
#  panel.grid=element_blank(),
#  panel.background=element_blank(),
#  panel.border=element_blank(),
#  axis.line=element_line(colour="black"),
#  plot.title=element_text(size=rel(1.2)))
#theme_set(old)
#----figure bar to KEGG ----
if(dim(KEGG[which(KEGG$P.Value < 0.05),])[1] < 30 | dim(KEGG[which(KEGG$P.Value < 0.05),])[1] > 50 ){
  KEGG <- KEGG[1:30,c("Name","P.Value")]  
}else{
  KEGG <- KEGG[which(KEGG$P.Value < 0.05),]
}
KEGG$Name <-as.character(KEGG$Name)
KEGG$Name <- as.factor(KEGG$Name)
KEGG$P.Value <- -log10(KEGG$P.Value)
p <- ggplot(KEGG,aes(x=reorder(Name,P.Value),
                         y=P.Value,
                         fill=P.Value))+
  geom_bar(stat="identity",position="dodge",width=0.6,colour="black")+
  coord_flip()+
  scale_fill_gradientn(colours = heat.colors(10)[10:1])+
  labs(title=name,fill="-log10(qvalue)")+xlab("")+ylab("-log10(qvalue)") + theme_classic()
#----save figure----
ggsave(filename=paste(outdir,"/",name,".KEGG.pdf",sep=""),plot=p,width = 8,height = 6)
ggsave(filename=paste(outdir,"/",name,".KEGG.png",sep=""),type="cairo-png",plot=p,width = 8,height = 6)


#----figure bar to GO ----
if(dim(GO[which(GO$P.Value < 0.05),])[1] < 30 | dim(GO[which(GO$P.Value < 0.05),])[1] > 50){
  GO <- GO[1:30,c("Name","P.Value")]  
}else{
  GO <- GO[which(GO$P.Value < 0.05),]
}
GO$Name <-as.character(GO$Name)
GO$Name <- as.factor(GO$Name)
GO$P.Value <- -log10(GO$P.Value)
p <- ggplot(GO,aes(x=reorder(Name,P.Value),
                     y=P.Value,
                     fill=P.Value))+
  geom_bar(stat="identity",position="dodge",width=0.6,colour="black")+
  coord_flip()+
  scale_fill_gradientn(colours = heat.colors(10)[10:1])+
  labs(title=name,fill="-log10(Corrected P-Value)")+xlab("")+ylab("-log10(Corrected P-Value)") + theme_classic()
#----save figure----
ggsave(filename=paste(outdir,"/",name,".GO.pdf",sep=""),plot=p,width = 8,height = 6)
ggsave(filename=paste(outdir,"/",name,".GO.png",sep=""),type="cairo-png",plot=p,width = 8,height = 6)


























args<-commandArgs(TRUE)
## for test
# name = 'haha'
# data.filepath = 'm1-12dpa_vs_m1-9dpa.down.KEGG.enrich.xls'
# enrich_type = 'down.KEGG'
# outdir = getwd()

name = args[1]
data.filepath = args[2]
enrich_type = args[3]
outdir = args[4]

library(ggplot2)

enrich_table <- read.delim(data.filepath,head=T,stringsAsFactors=F)

if (sum(grep("GO",enrich_type)) == 1) {
  names(enrich_table)[3] <- "qvalue"
  names(enrich_table)[6] <- "Name" 
} else {
  names(enrich_table)[1] <- "Name"
  names(enrich_table)[7] <- "qvalue"  
}

for(i in 1:dim(enrich_table)[1]){
  if(nchar(enrich_table$Name[i]) > 50){
    enrich_table$Name[i] <- substr(enrich_table$Name[i],1,50)
  }
}

#----figure bar to KEGG ----
enrich_table<- enrich_table[which(enrich_table$qvalue < 1),]

if(dim(enrich_table)[1] > 30 ){
  enrich_table <- enrich_table[1:30,c("Name","qvalue")]  
}
enrich_table$Name <-as.character(enrich_table$Name)
enrich_table$Name <- as.factor(enrich_table$Name)
enrich_table$qvalue <- -log10(enrich_table$qvalue)

enrich_num = dim(enrich_table)[1]
enrich_bin_with = ceiling(enrich_num/5)*0.1


if (enrich_num == 1) {
  p <- ggplot(enrich_table,aes(x=reorder(Name,qvalue),
                               y=qvalue))+
    geom_bar(stat="identity",position="dodge",width=enrich_bin_with,colour="black",fill="red")+
    coord_flip()+theme_classic()
} else {
  p <- ggplot(enrich_table,aes(x=reorder(Name,qvalue),
                               y=qvalue,
                               fill=qvalue))+
    geom_bar(stat="identity",position="dodge",width=enrich_bin_with,colour="black")+
    coord_flip()+scale_fill_gradientn(colours = heat.colors(10)[10:1])+theme_classic()
}
  
if (sum(grep("GO",enrich_type))) {
  p<- p+labs(title=name,fill="-log10(qvalue)")+xlab("")+ylab("-log10(qvalue)")
} else {
  p<- p+labs(title=name,fill="-log10(Corrected P-Value)")+xlab("")+ylab("-log10(Corrected P-Value)")
}

#----save figure----
ggsave(filename=paste(outdir,"/",name,".",enrich_type,".bar.pdf",sep=""),plot=p,width = 8,height = 6)
ggsave(filename=paste(outdir,"/",name,".",enrich_type,".bar.png",sep=""),type="cairo-png",plot=p,width = 8,height = 6)

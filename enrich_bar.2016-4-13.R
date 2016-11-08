args<-commandArgs(TRUE)
## for test
# name = 'haha'
# data.filepath = 'T3_dpi_col_vs_T6_dpi_col.all.KEGG.enrich.xls'
# enrich_type = 'KEGG'
# outdir = getwd()
# setwd("C:\\Users\\Administrator\\Desktop\\R\\data\\working")
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
enrich_table$qvalue <- -log10(enrich_table$qvalue)
enrich_table$Color <- "NA"
gradients <- function(data,locations,colors){
  n2 <- locations[2]
  n1 <- locations[1]
  if(length(colors) >= (n2-n1+1)){
    data[n1:n2] <- colors[1:(n2-n1+1)]
  }else{
    reps <-floor((n2-n1+1)/length(colors)) 
    if(n2-n1+1 > reps*length(colors)){
      data[n1:(n1+reps*length(colors)-1)] <- rep(colors,each=reps)
      data[(n1+reps*length(colors)):n2] <- colors[length(colors)]
    }else{
      data[n1:(n1+reps*length(colors)-1)] <- rep(colors,each=reps)
    }
  }
  return(data)
}
#设置颜色指标
up <- -log10(.001)
down <- -log10(.05)
up_l <- which(enrich_table$qvalue < up)[1]
down_l <- which(enrich_table$qvalue < down)[1]
#存在NA的情况
if(is.na(up_l) == T){  
  enrich_table$Color <- heat.colors(10)[1]  #q.value 全部大于 up
}else if(is.na(down_l) == T){
  if(up_l == 1){            #q.value 全部大于down 但全部小于up
    enrich_table$Color <- gradients(enrich_table$Color,c(1,dim(enrich_table)[1]),colors = heat.colors(15)[4:8])
  }else{        #q.value 全部大于down 但存在大于up
    enrich_table$Color <- gradients(enrich_table$Color,c(1,up_l-1),colors = heat.colors(15)[1])
    enrich_table$Color <- gradients(enrich_table$Color,c(up_l,dim(enrich_table)[1]),colors = heat.colors(15)[4:8])
  }
}else if(down_l == 1){           #q.value 全部小于 down
  enrich_table$Color <- gradients(enrich_table$Color,c(1,dim(enrich_table)[1]),colors = heat.colors(15)[10:15])
}else if(up_l == 1){
  if(is.na(down_l) == T){
    enrich_table$Color <- gradients(enrich_table$Color,c(1,dim(enrich_table)[1]),colors = heat.colors(15)[4:8])
  }else{
    enrich_table$Color <- gradients(enrich_table$Color,c(1,down_l-1),colors = heat.colors(15)[4:8])
    enrich_table$Color <- gradients(enrich_table$Color,c(down_l,dim(enrich_table)[1]),colors = heat.colors(15)[10:15])
  }
}else{                                #一般情况
  enrich_table$Color <- gradients(enrich_table$Color,c(1,up_l-1),colors = heat.colors(15)[1])
  enrich_table$Color <- gradients(enrich_table$Color,c(up_l,down_l-1),colors = heat.colors(15)[4:8])
  enrich_table$Color <- gradients(enrich_table$Color,c(down_l,dim(enrich_table)[1]),colors = heat.colors(15)[10:15])
}



#----figure bar to KEGG ----
enrich_table<- enrich_table[which(enrich_table$qvalue > 0),]

if(dim(enrich_table)[1] > 30 ){
  enrich_table <- enrich_table[1:30,c("Name","qvalue","Color")]  
}
# enrich_table$Name <-as.character(enrich_table$Name)
# enrich_table$Name <- as.factor(enrich_table$Name)

enrich_num = dim(enrich_table)[1]
enrich_bin_with = ceiling(enrich_num/5)*0.1
ymax <- ceiling(max(enrich_table$qvalue))
enrich_table_back <- enrich_table
enrich_table$qvalue  <- ifelse(enrich_table$qvalue > 3,
                               3+(enrich_table$qvalue - 3)/ymax,
                               enrich_table$qvalue)
enrich_table$Name<-factor(enrich_table$Name,levels = enrich_table$Name[dim(enrich_table)[1]:1])

heat_colors <- heat.colors(15)
names(heat_colors) <- heat.colors(15)

p <- ggplot(enrich_table,aes(x=Name,
                     y=qvalue,fill=Color))+
  geom_bar(stat="identity",position="dodge",width=enrich_bin_with,colour = "black")+
  scale_fill_manual(values = heat_colors)+
  coord_flip()+guides(fill = F)+theme_classic()

if (ymax <= 3){
  p<-p+scale_y_continuous(breaks = c(0:ymax),labels = c(0:ymax),limits = c(0,ymax))
} else {
  p<-p+scale_y_continuous(breaks = c(0:4),labels = c(0:3,ymax),limits = c(0,4))
} 

# p
# if (enrich_num == 1) {
#   p <- ggplot(enrich_table,aes(x=reorder(Name,qvalue),
#                                y=qvalue,fill=Color))+
#     geom_bar(stat="identity",position="dodge",width=enrich_bin_with,colour="black",fill="red")+
#     coord_flip()+theme_classic()
# } else {
#   p <- ggplot(enrich_table,aes(x=reorder(Name,qvalue),
#                                y=qvalue,
#                                fill=qvalue))+
#     geom_bar(stat="identity",position="dodge",width=enrich_bin_with,colour="black")+
#     coord_flip()+scale_fill_gradientn(colours = heat.colors(10)[10:1])+theme_classic()
# }

if (sum(grep("GO",enrich_type))) {
  p<- p+labs(title=name,fill="-log10(qvalue)")+xlab("")+ylab("-log10(qvalue)")
} else {
  p<- p+labs(title=name,fill="-log10(Corrected P-Value)")+xlab("")+ylab("-log10(Corrected P-Value)")
}

#----save figure----
ggsave(filename=paste(outdir,"/",name,".",enrich_type,".bar.pdf",sep=""),plot=p,width = 8,height = 6)
ggsave(filename=paste(outdir,"/",name,".",enrich_type,".bar.png",sep=""),type="cairo-png",plot=p,width = 8,height = 6)

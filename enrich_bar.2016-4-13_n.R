args<-commandArgs(TRUE)

## for test
# name = 'haha'
# data.filepath = 'T12_dpi_col_vs_T1_dpi_col.all.KEGG.enrich.xls'
# enrich_type = 'KEGG'
# outdir = getwd()

name = args[1]
data.filepath = args[2]
enrich_type = args[3]
outdir = args[4]

library(ggplot2)
library(ggthemes)

enrich_table <- read.delim(data.filepath,head=T,stringsAsFactors=F)

if (sum(grep("GO",enrich_type)) == 1) {
  names(enrich_table)[3] <- "qvalue"
  names(enrich_table)[6] <- "Name" 
} else {
  names(enrich_table)[1] <- "Name"
  names(enrich_table)[7] <- "qvalue"  
}

for(i in 1:dim(enrich_table)[1]){
  if(nchar(enrich_table$Name[i]) > 70){
    enrich_table$Name[i] <- paste(substr(enrich_table$Name[i],1,50),
substr(enrich_table$Name[i],nchar(enrich_table$Name[i])-10,nchar(enrich_table$Name[i])),sep = "...")
  }
}

enrich_table$qvalue[enrich_table$qvalue == 0] <- 1e-200
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

qvalues <- round(enrich_table$qvalue,1)
revse_heat_color <- heat.colors(31)[31:1]
names(revse_heat_color) <- seq(0,3,by = .1)
enrich_table$Color <- revse_heat_color[match(as.character(qvalues),names(revse_heat_color),nomatch = NA)]
enrich_table$Color[is.na(enrich_table$Color)] <- heat.colors(1)


#----figure bar to KEGG ----
enrich_table<- enrich_table[which(enrich_table$qvalue > 0),]

if (dim(enrich_table)[1] > 0){
  if(dim(enrich_table)[1] > 30 ){
    enrich_table <- enrich_table[1:30,c("Name","qvalue","Color")]  
  }

  enrich_num = dim(enrich_table)[1]
  enrich_bin_with = ceiling(enrich_num/5)*0.1
  ymax <- ceiling(max(enrich_table$qvalue))
  enrich_table_back <- enrich_table
  enrich_table$qvalue  <- ifelse(enrich_table$qvalue > 3,
                                 3+(enrich_table$qvalue - 3)/ymax,
                                 enrich_table$qvalue)
  enrich_table$Name<-factor(enrich_table$Name,levels = enrich_table$Name[dim(enrich_table)[1]:1])

  heat_colors <- enrich_table$Color
  names(heat_colors) <- enrich_table$Color

  p <- ggplot(enrich_table,aes(x=Name,
                               y=qvalue,fill=Color))+
    geom_bar(stat="identity",position="dodge",width=enrich_bin_with,colour = "black")+
    scale_fill_manual(values = heat_colors)+
    coord_flip()+guides(fill = F)+theme_few()

  if (ymax <= 3){
    p<-p+scale_y_continuous(breaks = c(0:ymax),labels = c(0:ymax),limits = c(0,ymax))
  } else {
    p<-p+scale_y_continuous(breaks = c(0:4),labels = c(0:3,ymax),limits = c(0,4))
  } 

  if (sum(grep("GO",enrich_type))) {
    p<- p+labs(title=name,fill="-log10(qvalue)")+xlab("")+ylab("-log10(qvalue)")
  } else {
    p<- p+labs(title=name,fill="-log10(Corrected P-Value)")+xlab("")+ylab("-log10(Corrected P-Value)")
  }

  #----save figure----
  ggsave(filename=paste(outdir,"/",name,".",enrich_type,".bar.pdf",sep=""),plot=p,width = 8,height = 6)
  ggsave(filename=paste(outdir,"/",name,".",enrich_type,".bar.png",sep=""),type="cairo-png",plot=p,width = 8,height = 6)
} else {
  print("no term/pathway qvalue smaller than 1")
}

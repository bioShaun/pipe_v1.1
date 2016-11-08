#2016-7-5 for enrichment barplot
#changed on 2016-10-14
library(ggplot2)
args<-commandArgs(TRUE)
options(stringsAsFactors = F)

enrich_type <- args[1]
name <- args[2]
go_enrich1 <- args[3]
go_enrich2 <- args[4]
diff_list1 <- args[5]
diff_list2 <- args[6]
go_file <- args[7]
out_dir <- args[8]

# enrich_type = "GO"
# name <- "A-a_vs_A-b"
name1  <- unlist(strsplit(name,"_vs_"))[1]
name2 <- unlist(strsplit(name,"_vs_"))[2]
GO.A_a <- read.delim(go_enrich1)
GO.A_b <- read.delim(go_enrich2)
GO.subsetA <- read.delim(diff_list1, header = F)
GO.subsetB <- read.delim(diff_list1, header = F )
GO.all <- read.delim(go_file, sep = ",")

# GO.A_a <- read.delim("A-a_vs_A-b.A-a-UP.GO.enrich.xls")
# GO.A_b <- read.delim("A-a_vs_A-b.A-b-UP.GO.enrich.xls")
# GO.subsetA <- read.delim("A-a-UP.subset.txt")
# GO.subsetB <- read.delim("A-b-UP.subset.txt")
# GO.all <- read.delim("tmp.Triticum_aestivum.go.txt",sep = ",")
term_count <- length(unique(GO.all[,c(1)]))
GO.A_count <- length(unique(GO.subsetA$V1))
GO.B_count <- length(unique(GO.subsetB$V1))

#cut overlong name
drop_overlong_name <- function(data,n=15){
  for(i in 1:dim(data)[1]){
    if(nchar(data$term[i]) > n){
      data$term[i] <- paste0(substr(data$term[i],1,n),'...')
    }
  }
  return(data)
}

#select data
data_select <- function(data){
  if(length(which(data$pvalue <= 0.05)) > 20){
    data <- data[which(data$pvalue <= 0.05)[1:20],]
  }else{
    data <- data[which(data$pvalue <= 0.05),]
  }
  return(data)
}

#prepare data
enrich_data <- function(data){
  data <- data[,c("over_represented_pvalue","numDEInCat","numInCat","term","ontology")]
  names(data)[1] <- 'pvalue'
  data <- data_select(drop_overlong_name(data))
  return(data)
}
data_setA <- enrich_data(GO.A_a)
data_setB <- enrich_data(GO.A_b)

#barpolt
library(ggthemes)
theme_set(theme_calc()+theme(
  axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5)
))
data_setA$color <- 'blue'
data_setB$color <- 'red'
label_a <- paste("No.",name1," up-regulated genes",sep = " ")
label_b <- paste("No.",name2," up-regulated genes",sep = " ")

all_data <- rbind(data_setA,data_setB) 
all_data <- all_data[,c("term","numDEInCat","numInCat","ontology","color")]
all_data$term <- factor(all_data$term,levels = all_data$term)

all_data$expected <- all_data$numInCat / term_count * c(GO.A_count,GO.B_count)
#-- add black color
all_data[dim(all_data)[1]+1,] <- all_data[dim(all_data)[1],]
all_data[dim(all_data)[1],c(2,6)] <- 0
all_data[dim(all_data)[1],5] <- 'black'
#--
col <- all_data$color
names(col) <- all_data$color

barplot <- ggplot(data = all_data,aes(x=term,y=numDEInCat,fill=color))+
  scale_fill_manual(values = col,labels = c('Expected no. of genes',label_a,label_b))+geom_bar(stat = 'identity', width = 0.75)+
  geom_bar(aes(x=term,y=expected),fill = 'black',stat = 'identity', width = 0.75)+
  facet_grid(.~ontology,scales = 'free_x',space = 'free')+guides(fill = guide_legend(title = ""))+
  xlab("") + ylab("Number of genes")

ggsave(filename = paste(out_dir, "/", enrich_type, "_", name, '_barplot.pdf',sep=""),plot = barplot,width = 0.6*(dim(all_data)[1] - 1),height = 8)
ggsave(filename = paste(out_dir, "/", enrich_type, "_", name, '_barplot.png',sep=""),type = 'cairo-png',plot = barplot,width = 0.6*(dim(all_data)[1] - 1),height = 8)

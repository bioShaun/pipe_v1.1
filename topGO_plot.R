library(topGO)
library(Cairo)
options(bitmapType='cairo')

args<-commandArgs(T)
gene_go_mape <- args[1]
diff_gene <- args[2]
enrich_result <- args[3]
name <- args[4]
out_dir <- args[5]

geneID2GO <- readMappings(file = gene_go_mape)
geneNames <- names(geneID2GO)
diff_gene <- read.table(file = diff_gene,header = F)
enrich_result <- read.table(file = enrich_result,header = T, sep = '\t', row.names = 1)
geneList <- factor(as.integer(geneNames %in% diff_gene$V1))
names(geneList) <- geneNames

go_catogary_vector <- c('MF','CC','BP')

bak_enrich_result <- enrich_result
enrich_result <- enrich_result[which(enrich_result$numInCat >= 5),]

for (i in 1:length(go_catogary_vector)) {
    go_catogary <- go_catogary_vector[i]
    each_enrich_result <- enrich_result[which(enrich_result$ontology == go_catogary),]
    each_go_qvalue <- each_enrich_result$qvalue
    each_go_qvalue[which(each_go_qvalue == 0)] <- 1e-100
    names(each_go_qvalue)<-row.names(each_enrich_result)
    if (dim(each_enrich_result)[1] < 2){
      out_info <- paste('Too little gene annotated to ',go_catogary,sep="")
      print(out_info)
    } else {
      GOdata <- new("topGOdata", ontology = go_catogary, allGenes = geneList,
        annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize = 5)  
      if (dim(each_enrich_result)[1] <= 10) {
        pdf(file = paste(out_dir,'/',name,'.',go_catogary,'.GO.DAG.pdf',sep = ''), 
            width = 8, height = 8)
        showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 1, useInfo = 'all')
        dev.off()
        Cairo(file=paste(out_dir,'/',name,'.',go_catogary,'.GO.DAG.png',sep = ''), 
          type="png",
          units="in", 
          width=8, 
          height=8, 
          pointsize=12, 
          dpi=300,
          bg = "white")
        showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 1, useInfo = 'all')
        dev.off()   
      } else { 
        pdf(file = paste(out_dir,'/',name,'.',go_catogary,'.GO.DAG.pdf',sep = ''), 
            width = 8, height = 8)
        showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 5, useInfo = 'all')
        dev.off()        
        Cairo(file=paste(out_dir,'/',name,'.',go_catogary,'.GO.DAG.png',sep = ''), 
          type="png",
          units="in", 
          width=8, 
          height=8, 
          pointsize=12, 
          dpi=300,
          bg = "white")
        showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 5, useInfo = 'all')
        dev.off()                 
      }
    ## plot pdf
    }
    # save(list=ls(all=TRUE), file='test.Rdata') 
}

# mf_enrich_result <- enrich_result[which(enrich_result$ontology == 'MF'),]
# mf_go_qvalue <- mf_enrich_result$qvalue
# names(mf_go_qvalue)<-row.names(mf_enrich_result)
# GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
#               annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize = 5)
# pdf(file = 'mf_go.pdf')
# showSigOfNodes(GOdata, mf_go_qvalue, firstSigNodes = 5, useInfo = 'all')
# dev.off()



# bp_enrich_result <- enrich_result[which(enrich_result$ontology == 'BP'),]
# bp_go_qvalue <- bp_enrich_result$qvalue
# names(bp_go_qvalue)<-row.names(bp_enrich_result)
# GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
#               annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize = 5)
# pdf(file = 'bp_go.pdf')
# showSigOfNodes(GOdata, bp_go_qvalue, firstSigNodes = 5, useInfo = 'all')
# dev.off()

# cc_enrich_result <- enrich_result[which(enrich_result$ontology == 'CC'),]
# cc_go_qvalue <- cc_enrich_result$qvalue
# names(cc_go_qvalue)<-row.names(cc_enrich_result)
# GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList,
#               annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize = 5)
# pdf(file = 'cc_go.pdf')
# showSigOfNodes(GOdata, cc_go_qvalue, firstSigNodes = 5, useInfo = 'all')
# dev.off()

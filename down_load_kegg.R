library(KEGGREST)
library(png)

args<-commandArgs(T)
spe <- args[1]
all_data_dir <- '/home/public/database/kegg/pathway/'
database_dir <- paste(all_data_dir,spe,sep = '/')
if (! dir.exists(database_dir)) { dir.create(database_dir) }

spe_pathinfo<-keggList("pathway", spe)
spe_pathid_info <- names(spe_pathinfo)
spe_pathid <- strsplit(spe_pathid_info,'path:')

for (i in 1:length(spe_pathid)) {
  pathway_id <-spe_pathid[[i]][2]
  xml <-keggGet(pathway_id, "kgml")
  writeChar(xml,con = paste(database_dir,'/',pathway_id,'.xml',sep=''),eos = NULL)
  png <- keggGet(pathway_id, "image")
  writePNG(png, paste(database_dir,'/',pathway_id,'.png',sep=''))
}

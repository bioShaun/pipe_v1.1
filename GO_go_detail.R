library(GO.db)
args<-commandArgs(T)
go_input <- args[1]
go_output <- args[2]

go_file <- read.table(go_input, header = T, sep = ',', stringsAsFactors= F)  
go_file <- subset(go_file, go_id != "")
go_list <- unique(go_file$go_id)
go_out <- as.data.frame(go_list)
go_term <- Term(go_list)
go_ontology <- Ontology(go_list)
go_out$go_ontology <- go_ontology  
go_out$go_term <- go_term     
write.table(go_out, file = go_output, quote = F,  sep = "\t", row.names = F)   

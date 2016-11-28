# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA);
args<-commandArgs(TRUE)
exp_data <- args[1]
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
myData = read.table(exp_data, header = T, row.names = 1, sep= '\t');
datExpr0 = as.data.frame(t(myData));
datExpr = datExpr0
save(datExpr, file = "dataInput.RData")
enableWGCNAThreads()
#lnames = load(file = "dataInput.RData");
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)
# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
pdf(file = "module_colors.pdf", width = 12, height = 9);
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "networkConstruction-auto.RData")

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
#modules = c("brown", "red");
# Select module probes
probes = names(datExpr)
#inModule = is.finite(match(moduleColors, modules));
#modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
#modTOM = TOM[inModule, inModule];
#dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(TOM,
  edgeFile = "CytoscapeInput-edges.txt",
  nodeFile = "CytoscapeInput-nodes.txt",
  weighted = TRUE,
  threshold = 0.1,
  nodeNames = probes,
  nodeAttr = moduleColors);

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
#sizeGrWindow(9,9)
pdf(file = "Network_heatmap.pdf", width = 12, height = 9);
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
#weight = as.data.frame(datTraits$weight_g);
#names(weight) = "weight"
# Add the weight to existing module eigengenes
#MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
#sizeGrWindow(5,7.5);
pdf(file = "plotEigengeneNetworks.pdf", width = 12, height = 9);
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()

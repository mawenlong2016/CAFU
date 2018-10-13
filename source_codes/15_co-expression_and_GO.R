#### Command to run tool:
# Rscript /home/malab10/software/galaxy/tools/myTools/WGCNA_and_topGO.R --input1 $input_1 --input2 $input_2 --input3 $input_3 --output1 $output_1 --output2 $output_2 --output3 $output_3 --output4 $output_4 --output5 $output_5
# Setup R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# We need to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Import required libraries
library("getopt")
options(stringAsfactors = F, useFancyQuotes = F)

# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Get options using the spec as defined by the enclosed list
# Read the options from the default: commandArgs(TRUE)
option_specification = matrix(c('input1', 'i1', 2, 'character',
                                'input2', 'i2', 2, 'character',
                                'input3', 'i3', 2, 'character',
                                'input4', 'i4', 2, 'character',
                                'input5', 'i5', 2, 'character',
                                'output1', 'o1', 2, 'character',
                                'output2', 'o2', 2, 'character',
                                'output3', 'o3', 2, 'character',
                                'output4', 'o4', 2, 'character',
                                'output5', 'o5', 2, 'character'
                               ),
                              byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

options(warn=-1)
library(WGCNA)
library(raster)
library(rsgcc)

## FPKM expression 

FPKMMat <- as.matrix(read.table(file = options$input1, header = T, sep = "\t", quote = ""))
class(FPKMMat) <- "numeric"

CVvalue <- options$input2
minExpSampleNum <- options$input3
coefficient <- options$input4
thread <- options$input5

FPKMMat[is.na(FPKMMat)] <- 0

numzero <- function(x){
  return(apply(x, 1, function(y){sum(y >= 1)}))
}
cvnum <- apply(FPKMMat,1,cv)/100

FPKMMat <- FPKMMat[cvnum>CVvalue,]
zeronum <- numzero(FPKMMat)
table(zeronum)
FPKMMat <- FPKMMat[which(zeronum >= minExpSampleNum),]
rm(cvnum, zeronum, numzero)

powers=c(c(1:10),seq(from=12,to=30,by=2))

sft=pickSoftThreshold(t(FPKMMat),powerVector=powers)

idx <- which(sft$fitIndices[ , 2] >= 0.9 & sft$fitIndices[ , 5] <= 200)
softPower <- idx[1]

class(FPKMMat) <- "numeric"
similarity = cor.matrix(FPKMMat,
           cpus = 10,
           cormethod = coefficient,
           style = c("all.pairs"),
           var1.id = NA,
           var2.id = NA,
           pernum = 0,
           sigmethod = c("two.sided"),
           output = c("matrix"))
#similarity = cor(t(FPKMMat), method = "spearman")
adjacency <- (0.5 + 0.5*similarity$corMatrix)^softPower
rm(similarity, softPower)
collectGarbage()
distTOM <- TOMdist(adjMat = adjacency, TOMType = "signed", TOMDenom = "min", verbose = 1, indent = 0)
collectGarbage()

# Hierarchical clustering:
hierTOM <- hclust(as.dist(distTOM),method="average")
cutHeight <- 0.95
minModuleSize <- 30
dynamicMods = cutreeDynamic(dendro = hierTOM, distM = distTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize, cutHeight = cutHeight);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
dynamicMods <- mergeCloseModules(t(FPKMMat), dynamicColors, cutHeight = 0.15)
dynamicColors = dynamicMods$colors
table(dynamicColors)
length(table(dynamicColors))

# Plot the dendrograms and module colors
pdf(file= options$output1, wi=9, he=5)
layout(matrix(c(1:2), 2,1), heights = c(0.8, 0.2));
plotDendroAndColors(hierTOM,
                    dynamicColors,
                    dendroLabels=FALSE,
                    setLayout = FALSE, abHeight = cutHeight,
                    marAll = c(0.3, 4.5, 2, 0.2))
dev.off();

## The results;
#for(i in seq(0.7,0.3,-0.1)){
i <- 0.55
cyt = exportNetworkToCytoscape(adjacency,
                               edgeFile=options$output2, nodeFile=options$output3,
                               weighted = TRUE, threshold = i, nodeNames=rownames(FPKMMat),
                               nodeAttr = dynamicColors)



####  The network analysis  ####
## Filtered Nodes and Edges;
Nodes <- read.table(options$output3, header = T, sep = "\t",stringsAsFactors = F)
Edges <- read.table(options$output2, header = T, sep = "\t",stringsAsFactors = F)


## ME
MEList = moduleEigengenes(t(FPKMMat[Nodes[,1],]), colors = Nodes[,3])
MEs = MEList$eigengenes
colnames(MEs) <- gsub("ME","",colnames(MEs))

## the hub genes   ## 0.55
novel <-  Nodes[,1][grepl("Contig", Nodes[,1])]
hubTransID <- matrix("NA", 1, 2)
for(i in unique(Nodes[,3])){
  curHubID <- Nodes[which(Nodes[,3]%in%i),1]
  curHubID <- FPKMMat[curHubID,]
  curHubID <- cor(t(curHubID),MEs[,i])
  curHubID <- head(sort(curHubID[,1], decreasing = T), n = floor(0.15 * length(which(Nodes[,3]%in%i)))) #floor(0.1 * length(which(Nodes[,3]%in%i))))
  curHubTrans <- matrix("NA", length(names(curHubID)), 2)
  curHubTrans[ , 1] <- names(curHubID)
  curHubTrans[ , 2] <- i
  hubTransID <- as.matrix(rbind(hubTransID, curHubTrans))
  cat(i, "\n")
}
hubTransID <- hubTransID[-1, ]
write.table(hubTransID, file = options$output4, quote = F, sep="\t", col.names = F, row.names = F)


library(ggplot2)
runTopGO <- function(geneID, statistic = "fisher", algorithm = "classic", topNodes = 100, plot = TRUE){
  
  if(!require(biomaRt)){
    source("https://www.bioconductor.org/biocLite.R")
    biocLite("biomaRt")
    n}
  
  if(!require(topGO)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("topGO")
  }
  
  mart <- biomaRt::useMart(biomart = "plants_mart", dataset = "zmays_eg_gene", host = 'plants.ensembl.org')  #athaliana_eg_gene zmays_eg_gene
  GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
  GTOGO <- GTOGO[GTOGO$go_id != '', ]
  geneID2GO <- by(GTOGO$go_id, GTOGO$ensembl_gene_id, function(x) as.character(x))
  
  all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
  int.genes <- geneID
  int.genes <- intersect(int.genes, names(geneID2GO))
  int.genes <- factor(as.integer(all.genes %in% int.genes))
  names(int.genes) = all.genes
  
  go.obj.BP <- new("topGOdata", ontology='BP'
                   , allGenes = int.genes
                   , annot = annFUN.gene2GO
                   , gene2GO = geneID2GO
  )
  
  go.obj.MF <- new("topGOdata", ontology='MF'
                   , allGenes = int.genes
                   , annot = annFUN.gene2GO
                   , gene2GO = geneID2GO
  )
  
  go.obj.CC <- new("topGOdata", ontology='CC'
                   , allGenes = int.genes
                   , annot = annFUN.gene2GO
                   , gene2GO = geneID2GO
  )
  
  ##########retrieve the gene list related to a GO ID######################
  allGO.BP <- genesInTerm(object = go.obj.BP)
  allGO.MF <- genesInTerm(object = go.obj.MF)
  allGO.CC <- genesInTerm(object = go.obj.CC)
  
  #########retrive the significant GO terms
  results.BP <- runTest(go.obj.BP, algorithm = algorithm, statistic = statistic)
  results.tab.BP <- GenTable(object = go.obj.BP, elimFisher = results.BP, topNodes = topNodes)
  #results.tab.BP[which(results.tab.BP$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  
  results.MF <- runTest(go.obj.MF, algorithm = algorithm, statistic = statistic)
  results.tab.MF <- GenTable(object = go.obj.MF, elimFisher = results.MF, topNodes = topNodes)
  #results.tab.MF[which(results.tab.MF$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  
  results.CC <- runTest(go.obj.CC, algorithm = algorithm, statistic = statistic)
  results.tab.CC <- GenTable(object = go.obj.CC, elimFisher = results.CC, topNodes = topNodes)
  #results.tab.CC[which(results.tab.CC$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  
  if(plot){
    df <- data.frame(Category = c(rep("BP", topNodes), rep("CC", topNodes), rep("MF", topNodes)), 
                     x = c(results.tab.BP$Significant, results.tab.CC$Significant, 
                           results.tab.MF$Significant),
                     y = c(-log10(as.numeric(results.tab.BP$elimFisher)), 
                           -log10(as.numeric(results.tab.CC$elimFisher)), 
                           -log10(as.numeric(results.tab.MF$elimFisher))),
                     size = c(-log10(as.numeric(results.tab.BP$elimFisher)),
                              -log10(as.numeric(results.tab.CC$elimFisher)), 
                              -log10(as.numeric(results.tab.MF$elimFisher)))
    )
    
    kk <- ggplot(data = df, aes(x = x, y = y)) + 
      geom_point(aes(color = Category, size = size)) + 
      scale_size_continuous(range = c(2,10)) + 
      labs(x = "The number of significant genes", y = "The adjusted p-values for each GO term")
  }
  print(kk)
  results <- list(BP = results.tab.BP, CC = results.tab.CC, MF = results.tab.MF)
  results
}


GORes <- matrix("NA", 1, 8)
for( i in 1:length(unique(Nodes[ , 3]))){
  
  idx <- which(Nodes[ , 3] == unique(Nodes[ , 3])[i])
  idx <- intersect(idx, grep("transcript", Nodes[ , 1]))
  curKnownTrans <- as.matrix(Nodes[idx, 1])
  
  kk <- curKnownTrans
  curKnownGene <- curKnownTrans
  for(j in 1:nrow(curKnownTrans)){
    
    kk[j, 1] <- strsplit(curKnownTrans[j, 1], ":")[[1]][2]
    curKnownGene[j, 1] <- strsplit(kk[j, 1], "_")[[1]][1]
    
  }
  
  topGOResult <- runTopGO(geneID = curKnownGene)
  
  BPResult <- topGOResult$BP
  idx <- which(as.numeric(BPResult[ , 6]) <= 0.05)
  curGORes <- matrix("NA", length(idx), 8)
  curGORes[ , 1] <- unique(Nodes[ , 3])[i]
  curGORes[ , 2] <- "BP"
  curGORes[ , c(3:8)] <- as.matrix(BPResult[idx, ])
  GORes <- as.matrix(rbind(GORes, curGORes))
  
  CCResult <- topGOResult$CC
  idx <- which(as.numeric(CCResult[ , 6]) <= 0.05)
  curGORes <- matrix("NA", length(idx), 8)
  curGORes[ , 1] <- unique(Nodes[ , 3])[i]
  curGORes[ , 2] <- "CC"
  curGORes[ , c(3:8)] <- as.matrix(CCResult[idx, ])
  GORes <- as.matrix(rbind(GORes, curGORes))
  
  MFResult <- topGOResult$MF
  idx <- which(as.numeric(MFResult[ , 6]) <= 0.05)
  curGORes <- matrix("NA", length(idx), 8)
  curGORes[ , 1] <- unique(Nodes[ , 3])[i]
  curGORes[ , 2] <- "MF"
  curGORes[ , c(3:8)] <- as.matrix(MFResult[idx, ])
  GORes <- as.matrix(rbind(GORes, curGORes))

}
GORes <- GORes[-1, ]
write.table(GORes, options$output5, sep	= "\t", quote = F, col.names = F, row.names = F)

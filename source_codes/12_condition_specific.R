##Command to run tool:
# Rscript /mnt/xvdf/galaxy/tools/myTools/tissue_specifc.R --input1 $output --input2 $sampleInfo --input3 $TSVal--output1 $output_1 --output2 $output_2 ; 

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
                                'input3', 'o3', 2, 'character',
                                'input4', 'o4', 2, 'character',
                                'output1', 'o1', 2, 'character',
                                'output2', 'o2', 2, 'character'
                               ),
                              byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

# scripts
allFPKM <- as.matrix(read.table(file = options$input1, sep = "\t", quote = "", header = T, row.names = 1))
sampleInfo <- as.matrix(read.table(file = options$input2, sep = ";", quote = "", header = F))
TSVal <- options$input3
novelTransID <- as.matrix(read.table(file = options$input4, sep = "\t", quote = "", header = F))



colnames(allFPKM) <- sub("_exp.isoforms.results", "", colnames(allFPKM))
colnames(allFPKM) <- sub("/home/chenzhuod/galaxy/database/files/Expression/expression_evidence/", "", colnames(allFPKM))
rownames(allFPKM) <- sub("\"", "", rownames(allFPKM))
rownames(allFPKM) <- sub("\"", "", rownames(allFPKM))

##Tissue-specifc identity
novelExp <- allFPKM[novelTransID[ , 1], ]

sampleNum <- matrix("NA", nrow = nrow(sampleInfo), ncol = 4)
for( i in 1:nrow(sampleInfo)){
  
  curTissueID <- sampleInfo[i, 1]
  curSampleNum <- strsplit(sampleInfo[i, 2], ",")[[1]]
  
  sampleNum[i, 1] <- curTissueID
  sampleNum[i, 2] <- length(curSampleNum)
  
}
sampleID <- c()
for( i in 1:nrow(sampleInfo)){
  if( i == 1 ){
    sampleNum[i, 3] <- 1
    sampleNum[i, 4] <- sampleNum[i, 2]
  } else {
    
    sampleNum[i, 3] <- as.numeric(sampleNum[(i - 1), 4]) + 1
    sampleNum[i, 4] <- sum(as.numeric(sampleNum[c(1:i), 2]))
    
  }
  curSampleID <- strsplit(sampleInfo[i, 2], ",")[[1]]
  sampleID <- c(sampleID, curSampleID)
}
novelExp <- novelExp[ , sampleID]
medianExp <- matrix("NA", nrow = nrow(novelExp), nrow(sampleInfo))

for( i in 1:ncol(medianExp)){
  
  curExp <- novelExp[ , c(as.numeric(sampleNum[i, 3]):as.numeric(sampleNum[i, 4]))]
  class(curExp) <- "numeric"
  
  medianExp[ , i] <- apply(curExp, 1, median)	
  
}
rownames(medianExp) <- rownames(novelExp)
colnames(medianExp) <- sampleInfo[ , 1]

## tissue-specific score
scoreMat <- medianExp
for(i in 1:nrow(scoreMat)){
  for(j in 1:ncol(scoreMat)){
    scoreMat[i, j] <- 1 - (median(as.numeric(medianExp[i, -j])) / (as.numeric(medianExp[i, j])))
  }
  cat(i, j, "\n")
}

write.table(scoreMat, "/home/chenzhuod/CAFU_data/chenzhuod@gmail.com/CS/04scoreMat", sep = "\t", quote = F, col.names = T, row.names = T)


##filter tissue-specif transcripts
extractCSCandidate <- function(x){
  if(length(which(x > as.numeric(TSVal))) == 1){
    res <- "yes"
  } else {
    res <- "no"
  }
  res
}
class(scoreMat) <- "numeric"
tissueSpecificScore <- scoreMat[names(which(apply(scoreMat, 1, extractCSCandidate) == "yes")), ]

##
write.table(names(which(apply(scoreMat, 1, extractCSCandidate) == "yes")), "/home/chenzhuod/CAFU_data/chenzhuod@gmail.com/CS/05tissueSpecificID", sep = "\t", quote = F, col.names = T, row.names = T)
write.table(tissueSpecificScore, "/home/chenzhuod/CAFU_data/chenzhuod@gmail.com/CS/05tissueSpecificScore", sep = "\t", quote = F, col.names = T, row.names = T)


tissueSpecificID <- NULL
for(i in 1:ncol(tissueSpecificScore)){
  
  idx <- which(as.numeric(tissueSpecificScore[ , i]) > as.numeric(TSVal))
  curScore <- as.matrix(tissueSpecificScore[idx, ])
  curRes <- matrix("NA", nrow = nrow(curScore), 2)
  curRes[ , 1] <- rownames(curScore)
  curRes[ , 2] <- colnames(tissueSpecificScore)[i]
  
  tissueSpecificID <- as.matrix(rbind(tissueSpecificID, curRes))
  
  cat(colnames(tissueSpecificScore)[i], length(idx), "\n")
}
write.table(tissueSpecificID, file = options$output1, sep = "\t", quote = F, col.names = F, row.names = F)

##pheatmap
tisSpeTransExp <- allFPKM[tissueSpecificID[ , 1], ]

class(tisSpeTransExp) <- "numeric"
MeanVal <- apply(tisSpeTransExp, 1, mean)
SDVal <- apply(tisSpeTransExp, 1, sd)

Zscore <- (tisSpeTransExp - MeanVal)/SDVal
class(Zscore) <- "numeric"

library(pheatmap)

pdf(file = options$output2, height = 5, width = 10)
breaks <- unique(c(seq(min(Zscore), 3, 0.1), seq(3, max(Zscore), 0.5)))
pheatmap(Zscore, breaks = breaks, color = colorRampPalette(rev(c("red", "yellow")))(length(breaks)), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
dev.off()


library(ggplot2)
library(pheatmap)

allFPKM <- as.matrix(read.table(file = "/home/chenzhuod/CAFU_data/chenzhuod@gmail.com/CS/assemble_Exp", sep = "\t", quote = "", header = T, row.names = 1))
sampleInfo <- as.matrix(read.table(file = "/home/chenzhuod/CAFU_data/chenzhuod@gmail.com/CS/CS_RNA-Seq_Info.txt", sep = ";", quote = "", header = F))
TSVal <- 0.8
novelTransID <- as.matrix(read.table(file = "/home/chenzhuod/CAFU_data/chenzhuod@gmail.com/CS/635_novel_transcript_ID.txt", sep = "\t", quote = "", header = F))

colnames(allFPKM) <- sub("_exp.isoforms.results", "", colnames(allFPKM))
colnames(allFPKM) <- sub("/home/chenzhuod/galaxy/database/files/Expression/expression_evidence/", "", colnames(allFPKM))
rownames(allFPKM) <- sub("\"", "", rownames(allFPKM))
rownames(allFPKM) <- sub("\"", "", rownames(allFPKM))

## Tissue-specifc identity
novelExp <- allFPKM[novelTransID[ , 1], ]

medianExp <- matrix("NA", nrow = nrow(novelExp), nrow(sampleInfo))
rownames(medianExp) <- rownames(novelExp)
colnames(medianExp) <- sampleInfo[ , 1]

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

for( i in 1:ncol(medianExp)){
  
  curExp <- novelExp[ , c(as.numeric(sampleNum[i, 3]):as.numeric(sampleNum[i, 4]))]
  class(curExp) <- "numeric"
  
  medianExp[ , i] <- apply(curExp, 1, median)	
  
}

### tissue-specific score
scoreMat <- medianExp
for(i in 1:nrow(scoreMat)){
  for(j in 1:ncol(scoreMat)){
    scoreMat[i, j] <- 1 - (median(as.numeric(medianExp[i, -j])) / (as.numeric(medianExp[i, j])))
  }
  cat(i, j, "\n")
}

### filter tissue-specif transcripts
idx <- c()
for(i in 1:nrow(scoreMat)){
  if(length(which(as.numeric(scoreMat[i, ]) > TSVal)) == 1){
    idx <- c(idx, i)
  }
}
tissueSpecificScore <- scoreMat[idx, ]

tissueSpecificID <- NULL
for(i in 1:ncol(tissueSpecificScore)){
  
  idx <- which(as.numeric(tissueSpecificScore[ , i]) > TSVal)
  curScore <- as.matrix(tissueSpecificScore[idx, ])
  curRes <- matrix("NA", nrow = nrow(curScore), 2)
  curRes[ , 1] <- rownames(curScore)
  curRes[ , 2] <- colnames(tissueSpecificScore)[i]
  
  tissueSpecificID <- as.matrix(rbind(tissueSpecificID, curRes))
  
  cat(colnames(tissueSpecificScore)[i], length(idx), "\n")
}
write.table(tissueSpecificID, file = "/home/chenzhuod/CAFU_data/chenzhuod@gmail.com/CS/TS_ID", sep = "\t", quote = F, col.names = F, row.names = F)



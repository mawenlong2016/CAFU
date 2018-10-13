#### Command to run tool:
# Rscript /mnt/xvdf/galaxy/tools/myTools/integrate_gmapRes.R --output1 $gmapRes_integrated ;

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
option_specification = matrix(c('output1', 'o1', 2, 'character'
                               ),
                              byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

filterHitNum <- function(x){
  if((x[2] != "NA") & (x[3] == "NA")){
    return <- "Yes"
  } else {
    return <- "No"
  }
}

resDic <- "/home/chenzhuod/galaxy/database/files/Transcript_evidence/genome_level/gmapRes/"
resDir <- list.files(resDic)

allTransID <- as.matrix(read.table("/home/chenzhuod/galaxy/database/files/Transcript_evidence/genome_level/transcripts_ID", sep = "\t", quote = "", header = F))

gmapRes <- matrix("NA", nrow = nrow(allTransID), ncol = (length(resDir) * 2))
rownames(gmapRes) <- allTransID
colnames(gmapRes) <- rep(sub(".fa_gmapRes.psl", "", resDir), each = 2)

colID <- seq(1, ncol(gmapRes), 2)

for(i in 1:length(resDir)){
  
  curDir <- paste0(resDic, resDir[i])
  curGmapRes <- as.matrix(read.table(file = curDir, sep = "\t", quote = "", header = F))
  
  mapFrequency <- as.matrix(table(curGmapRes[ , 10]))
  maxFrequency <- max(mapFrequency[ , 1])
  
  coverage <- as.matrix((as.numeric(curGmapRes[ , 13]) - as.numeric(curGmapRes[ , 12])) / as.numeric(curGmapRes[ , 11]))
  identity <- as.matrix(as.numeric(curGmapRes[ , 1]) / (as.numeric(curGmapRes[ , 13]) - as.numeric(curGmapRes[ , 12])))
  curGmapRes <- as.matrix(cbind(coverage, identity, curGmapRes))
  
  curGmapInfo <- matrix("NA", nrow = nrow(curGmapRes), ncol = 3)
  curGmapInfo[ , c(1, 2)] <- curGmapRes[ , c(12, 3)]
  for(j in 1:nrow(curGmapRes)){
    curGmapInfo[j, 3] <- paste(as.numeric(curGmapRes[j, 1]), as.numeric(curGmapRes[j, 2]), as.numeric(curGmapRes[j, 3]), curGmapRes[j, 11], as.numeric(curGmapRes[j, 14]), as.numeric(curGmapRes[j, 15]), curGmapRes[j, 16], as.numeric(curGmapRes[j, 17]), as.numeric(curGmapRes[j, 18]), as.numeric(curGmapRes[j, 19]), sep = "_")
  }
  
  gmapInfo <- matrix("NA", nrow = nrow(mapFrequency), ncol = (maxFrequency + 1))
  gmapInfo[ , 1] <- rownames(mapFrequency)
  for(j in 1:nrow(gmapInfo)){
    if(mapFrequency[j, 1] == 1){
      kk <- t(curGmapInfo[which(curGmapInfo[ , 1] == gmapInfo[j, 1]), c(2, 3)])
      gmapInfo[j, 2] <- kk[ , 2] 
    }else{
      kk <- curGmapInfo[which(curGmapInfo[ , 1] == gmapInfo[j, 1]), c(2, 3)]
      tt <- kk[order(as.numeric(kk[ , 1]), decreasing = T), ]
      gmapInfo[j, c(2:(nrow(tt) + 1))] <- tt[c(1:nrow(tt)), 2]
    }
  }
  
  unmappedID <- setdiff(allTransID[ , 1], gmapInfo[ , 1])
  unmappedInfo <- matrix("NA", nrow = length(unmappedID), ncol = ncol(gmapInfo))
  unmappedInfo[ , 1] <- unmappedID
  
  allInfo <- as.matrix(rbind(gmapInfo, unmappedInfo))
  allInfo <- allInfo[order(allInfo[ , 1]), ]
  
  oneHitTransAnno <- allInfo[which(apply(allInfo, 1, filterHitNum) == "Yes"), ]
  
  for(j in 1:nrow(oneHitTransAnno)){
    
    gmapRes[oneHitTransAnno[j, 1], colID[i]] <- as.numeric(strsplit(oneHitTransAnno[j, 2], "_")[[1]][1])
    gmapRes[oneHitTransAnno[j, 1], (colID[i] + 1)] <- as.numeric(strsplit(oneHitTransAnno[j, 2], "_")[[1]][2])
    
  }

}

write.table(gmapRes, file = options$output1, sep = "\t", quote = F, col.names = T, row.names = T)

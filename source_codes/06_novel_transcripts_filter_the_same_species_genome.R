#### Command to run tool:
# Rscript /mnt/xvdf/galaxy/tools/myTools/novel_transcripts_filter_the_same_species_genome.R --output1 $same_similar_transciript_ID --output2 $novel_transcript_ID --input2 $Parameters.minCoverage --input3 $Parameters.minIdentity;

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
                                'output1', 'o1', 2, 'character',
                                'output2', 'o2', 2, 'character',
                                'output3', 'o3', 2, 'character'
                               ),
                              byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

# input gmap results and gene annotation
filteredTransID <- as.matrix(read.table(file = "/home/chenzhuod/galaxy/database/files/Transcript_evidence/genome_level/filteredTransID", sep = "\t", quote = "", header = F))
BestMultiIntronInfo <- as.matrix(read.table(file = "/home/chenzhuod/galaxy/database/files/Transcript_evidence/genome_level/BestMultiIntronInfo", sep = "\t", quote = "", header = F))
rownames(BestMultiIntronInfo) <- BestMultiIntronInfo[ , 1]
BestMultiIntronInfo <- BestMultiIntronInfo[ , -1]

BestMultiExon_knownTrans <- as.matrix(read.table(file = "/home/chenzhuod/galaxy/database/files/Transcript_evidence/genome_level/multiExonTrans_refTrans_bedtools_intersectRes", sep = "\t", quote = "", header = F))
BestMultiExonIntergenic <- BestMultiExon_knownTrans[which(as.numeric(BestMultiExon_knownTrans[ , 19]) == 0), ]
BestMultIExonIntergenicID <- as.matrix(BestMultiExonIntergenic[ , 9])
BestMultiExonCodingRegion <- BestMultiExon_knownTrans[which(as.numeric(BestMultiExon_knownTrans[ , 19]) != 0), ]

BestMutliOverlapID <- matrix("NA", nrow = nrow(BestMultiExonCodingRegion), ncol = 2)
BestMutliOverlapID[ , 1] <- BestMultiExonCodingRegion[ , 9]
kk <- matrix("NA", nrow = nrow(BestMultiExonCodingRegion), ncol = 1)
for(i in 1:nrow(BestMultiExonCodingRegion)){
  kk[i, 1] <- as.matrix(strsplit(BestMultiExonCodingRegion[i, 18], "ID=transcript:")[[1]][2])
  BestMutliOverlapID[i, 2] <- as.matrix(strsplit(kk[i, 1], ";")[[1]][1])
}

knwonIntonGFF <- as.matrix(read.table(file = "/home/chenzhuod/galaxy/database/files/Transcript_evidence/genome_level/all_intron_annotation", sep = "\t", quote = "", header = F))
knwonIntonGFFGeneID <- matrix("NA", nrow = nrow(knwonIntonGFF), ncol = 1)
knwonIntonGFFTransID <- matrix("NA", nrow = nrow(knwonIntonGFF), ncol = 1)
kk <- matrix("NA", nrow = nrow(knwonIntonGFF), ncol = 1)
for(i in 1:nrow(knwonIntonGFF)){
  kk[i, ] <- as.matrix(strsplit(knwonIntonGFF[i, 9], "ID=gene:")[[1]][2])
  knwonIntonGFFGeneID[i, ] <- as.matrix(strsplit(kk[i, 1], ";")[[1]][1])
  knwonIntonGFFTransID[i, ] <- as.matrix(strsplit(knwonIntonGFF[i, 9], "Parent=transcript:")[[1]][2])
  cat(i, "\n")
}
knwonIntonGFF <- cbind(knwonIntonGFF, knwonIntonGFFGeneID, knwonIntonGFFTransID)

KnownIntronInfo <- matrix("NA", nrow = nrow(as.matrix(unique(sort(knwonIntonGFF[ , 11])))), ncol = (max(table(knwonIntonGFF[ , 11])) + 2))
KnownIntronInfo[ , 1] <- as.matrix(unique(sort(knwonIntonGFF[ , 11])))
for(i in 1:nrow(KnownIntronInfo)){
  kk <- as.matrix(knwonIntonGFF[which(knwonIntonGFF[ , 11] == KnownIntronInfo[i, 1]), c(4, 5, 7)])
  if(ncol(kk) != 1){
    for(j in 1:nrow(kk)){
      KnownIntronInfo[i, 2] <- kk[1, 3]
      KnownIntronInfo[i, (j + 2)] <- as.matrix(paste(kk[j, 1], kk[j, 2], sep = ":"))
    }
  }else{
    KnownIntronInfo[i, 2] <- kk[3, 1]
    KnownIntronInfo[i, 3] <- as.matrix(paste(kk[1, 1], kk[2, 1], sep = ":"))
  }
}

curKnownIntronInfo <- KnownIntronInfo
curKnownIntronInfo <- sub(" ", "", KnownIntronInfo)
curKnownIntronInfo <- sub(": ", ":", curKnownIntronInfo)
rownames(curKnownIntronInfo) <- curKnownIntronInfo[ , 1]
curKnownIntronInfo <- curKnownIntronInfo[ , c(-1, -2)]

rownames(BestMutliOverlapID) <- BestMutliOverlapID[ , 2]
curBestMultiOverlapID <- matrix("NA", 1, ncol = ncol(BestMutliOverlapID))
for(i in 1:nrow(BestMutliOverlapID)){
  if(BestMutliOverlapID[i, 2] %in% rownames(curKnownIntronInfo)){
    curBestMultiOverlapID <- rbind(curBestMultiOverlapID, BestMutliOverlapID[i, ])
  }else{
    curBestMultiOverlapID <- curBestMultiOverlapID
  }
}
curBestMultiOverlapID <- curBestMultiOverlapID[-1, ]

BestMultiJudgeResult <- matrix("NA", nrow = nrow(curBestMultiOverlapID), 6)
rownames(BestMultiJudgeResult) <- curBestMultiOverlapID[ , 1]
colnames(BestMultiJudgeResult) <- c("assTransIntronNum", "mapRefTrans", "refTransIntronNum", "diffLessThan10Num", "diffEqual0Num", "")

for(i in 1:nrow(curBestMultiOverlapID)){
  ## assembled transcript intron number
  BestMultiJudgeResult[i, 1] <- length(BestMultiIntronInfo[curBestMultiOverlapID[i, 1], c(1:(ncol(BestMultiIntronInfo) - (length(which(BestMultiIntronInfo[curBestMultiOverlapID[i, 1], ] == "NA")))))])
  ## overlap (map) ref transcript
  BestMultiJudgeResult[i, 2] <- curBestMultiOverlapID[i, 2]
  ## ref transcript intron number
  BestMultiJudgeResult[i, 3] <- length(curKnownIntronInfo[curBestMultiOverlapID[i, 2], c(1:(ncol(curKnownIntronInfo) - (length(which(curKnownIntronInfo[curBestMultiOverlapID[i, 2], ] == "NA")))))])
  
  ## assemble transcript intron info
  NIntronInfo <- as.matrix(BestMultiIntronInfo[curBestMultiOverlapID[i, 1], c(1:(ncol(BestMultiIntronInfo) - (length(which(BestMultiIntronInfo[curBestMultiOverlapID[i, 1], ] == "NA")))))])
  ## ref transcript intron info
  KIntronInfo <- as.matrix(curKnownIntronInfo[curBestMultiOverlapID[i, 2], c(1:(ncol(curKnownIntronInfo) - (length(which(curKnownIntronInfo[curBestMultiOverlapID[i, 2], ] == "NA")))))])
  
  ## difference of intron between assemble and reference
  if(nrow(NIntronInfo) <= nrow(KIntronInfo)){
    IntronCompare <- matrix("NA", nrow = nrow(NIntronInfo), ncol = nrow(KIntronInfo))
    for(j in 1:nrow(NIntronInfo)){
      for(k in 1:nrow(KIntronInfo)){
        IntronCompare[j, k] <- as.matrix((abs(as.numeric(strsplit(NIntronInfo[j, 1], ":")[[1]][1]) - as.numeric(strsplit(KIntronInfo[k, 1], ":")[[1]][1]))) + (abs(as.numeric(strsplit(NIntronInfo[j, 1], ":")[[1]][2]) - as.numeric(strsplit(KIntronInfo[k, 1], ":")[[1]][2]))))
      }
    }
    ## intron difference feature stat
    BestMultiJudgeResult[i, 4] <- length(which(as.numeric(IntronCompare) < 10))
    BestMultiJudgeResult[i, 5] <- length(which(as.numeric(IntronCompare) == 0))
    kk <- matrix("NA", 1, ncol = nrow(NIntronInfo))
    for(j in 1:nrow(NIntronInfo)){
      kk[1, j] <- min(as.numeric(IntronCompare[j, ]))
    }
    ##
    BestMultiJudgeResult[i, 6] <- sum(as.numeric(kk[1, which(as.numeric(kk[1, ]) < 10)]))
  }
}

sameTransID <- unique(rownames(BestMultiJudgeResult[which(as.numeric(BestMultiJudgeResult[ , 1]) == as.numeric(BestMultiJudgeResult[ , 5])), ]))
diffTrans <- BestMultiJudgeResult[which(BestMultiJudgeResult[ , 1] != BestMultiJudgeResult[ , 5]), ]
idx <- c()
for( i in 1:nrow(diffTrans)){
  if(rownames(diffTrans)[i] %in% sameTransID){
    idx <- c(idx, i)
  }
}

if(length(idx) != 0){
  diffTrans <- diffTrans[-idx, ]
}
similarID <- unique(rownames(diffTrans[which(as.numeric(diffTrans[ , 1]) == as.numeric(diffTrans[ , 4])), ]))

same_similar_transciript_ID <- unique(union(sameTransID, similarID))
novelTransID <- setdiff(filteredTransID, same_similar_transciript_ID)


write.table(same_similar_transciript_ID, file = options$output1, sep = "\t", quote = F, col.names = F, row.names = F)
write.table(novelTransID, file = options$output2, sep = "\t", quote = F, col.names = F, row.names = F)

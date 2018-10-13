#### Command to run tool:
# Rscript /mnt/xvdf/galaxy/tools/myTools/get_same_species_best_alignment_filtered.R --input $sameSpeciesRefID ;

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
								'input3', 'i3', 2, 'character'
                               ),
                              byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

# input gmap results and gene annotation
filterHitNum <- function(x){
  if((x[2] != "NA") & (x[3] == "NA")){
    return <- "Yes"
  } else {
    return <- "No"
  }
}

filterMappedTrans <- function( x ){
  coverage <- as.numeric(strsplit(x[2], "_")[[1]][1])
  identity <- as.numeric(strsplit(x[2], "_")[[1]][2])
  if( coverage >= minCoverage & identity >= minIdentity ){
    return <- "Yes"
  } else {
    return <- "No"
  }
}

minCoverage <- options$input2
minIdentity <- options$input3


sameSpeciesGmapResDir <- paste0("/home/chenzhuod/galaxy/database/files/Transcript_evidence/genome_level/gmapRes/", options$input1, ".fa_gmapRes.psl")
curGmapRes <- as.matrix(read.table(sameSpeciesGmapResDir, sep = "\t", quote = "", header = F))
rownames(curGmapRes) <- curGmapRes[ , 10]
allTransID <- as.matrix(read.table("/home/chenzhuod/galaxy/database/files/Transcript_evidence/genome_level/transcripts_ID", sep = "\t", quote = "", header = F))

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
mappedTransID <- oneHitTransAnno[apply(oneHitTransAnno, 1, filterMappedTrans) == "Yes", 1]

mappedTransGmapRes <- curGmapRes[as.matrix(mappedTransID), ]
rownames(mappedTransGmapRes) <- mappedTransGmapRes[ , 12]
idx <- which(as.numeric(mappedTransGmapRes[ , 20]) == 1)
multiExonTransGmapRes <- mappedTransGmapRes[names(mappedTransGmapRes[-idx, 10]), ]

multiExonTransGmapResGFF3 <- matrix("NA", nrow = nrow(multiExonTransGmapRes), ncol = 9)

multiExonTransGmapResGFF3[ , 4] <- as.matrix(as.numeric(multiExonTransGmapRes[ , 18]) + 1)
multiExonTransGmapResGFF3[ , 5] <- as.matrix(as.numeric(multiExonTransGmapRes[ , 19]))
multiExonTransGmapResGFF3[ , c(1, 7, 9)] <- multiExonTransGmapRes[ , c(16, 11, 12)]
multiExonTransGmapResGFF3[ , 2] <- "B73"
multiExonTransGmapResGFF3[ , 3] <- "gmap"
multiExonTransGmapResGFF3[ , c(6, 8)] <- "."

BestMultiIntronInfo <- matrix("NA", nrow = nrow(multiExonTransGmapRes), ncol = max(as.numeric(multiExonTransGmapRes[ , 20])))
BestMultiIntronInfo[ , 1] <- multiExonTransGmapRes[ , 12]
for(i in 1:nrow(BestMultiIntronInfo)){
  startLocus <- as.matrix(strsplit(multiExonTransGmapRes[i, 23], ",")[[1]])
  length <- as.matrix(strsplit(multiExonTransGmapRes[i, 21], ",")[[1]])
  if(as.numeric(multiExonTransGmapRes[i, 20]) != 1){
    for(j in 1:(nrow(startLocus) - 1)){
      BestMultiIntronInfo[i, (j + 1)] <- as.matrix(paste((as.numeric(startLocus[j, 1]) + as.numeric(length[j, 1]) + 1), as.numeric(startLocus[(j + 1), 1]), sep = ":"))
    }
  }
}

write.table(BestMultiIntronInfo, file = "/home/chenzhuod/galaxy/database/files/Transcript_evidence/genome_level/BestMultiIntronInfo", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(multiExonTransGmapResGFF3, file = "/home/chenzhuod/galaxy/database/files/Transcript_evidence/genome_level/multiExonTransGmapRes.gff3", sep = "\t", quote = F, col.names = F, row.names = F)

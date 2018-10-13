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
option_specification = matrix(c('output1', 'o1', 2, 'character',
								'output2', 'o2', 2, 'character'
                               ),
                              byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);


resDic <- "/home/chenzhuod/galaxy/database/files/Expression/expression_evidence/coverage/"
resDir <- list.files(resDic)

for( i in 1:length(resDir)){
  
  curDir <- paste0(resDic, resDir[i], "/")
  sampleDir <- list.files(curDir)
  
  if( i == 1 ){
    
    readCovRes <- matrix("NA", nrow = length(resDir), ncol = length(sampleDir))
    rownames(readCovRes) <- list.files(resDic)
    colnames(readCovRes) <- list.files(curDir)
  }
  
  for( j in 1:length(sampleDir)){
    
    curSamDir <- paste0(curDir, sampleDir[j])
    curSamRes <- as.matrix(read.table(file = curSamDir, sep = "\t", quote = "", header = F))
    
    if( j == 1 ){
      
      curReadCovRes <- matrix("NA", nrow = nrow(curSamRes), ncol = length(sampleDir))
      curReadCovRes[ , j] <- curSamRes[ , 5]
      
    } else {
      
      curReadCovRes[ , j] <- curSamRes[ , 5]
      
    }
    
    idx <- which(as.numeric(curReadCovRes[ , j]) >= 5)
    readCovRes[i, j] <- length(idx)/nrow(curReadCovRes)
    
  }

  cat(i, nrow(curReadCovRes), "\n")
}

write.table(readCovRes, file = options$output1, sep = "\t", quote = F, col.names = T, row.names = T)

class(readCovRes) <- "numeric"
idx <- which(apply(readCovRes, 1, max) > 0.8)

confTransID <- rownames(readCovRes[idx, ])
write.table(confTransID, file = options$output2, sep = "\t", quote = F, col.names = F, row.names = F)
  

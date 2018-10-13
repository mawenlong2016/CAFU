#### Command to run tool:
# Rscript /mnt/xvdf/galaxy/tools/myTools/filter_transcript.R --input1 $gmapRes_integrated --input2 $minCoverage --input3 $minIdentity --output1 $filteredTransInfo --output2 $filteredTransID;

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
								'output1', 'o1', 2, 'character',
								'output2', 'o2', 2, 'character'
                               ),
                              byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

gmapRes <- as.matrix(read.table(options$input1, sep = "\t", quote = "", header = T))
minCoverage <- options$input2
minIdentity <- options$input3


filteredTransInfo <- matrix("NA", 1, 2)

colId <- seq(1, ncol(gmapRes), 2)
for(i in 1:length(colId)){
  
  idx <- which((as.numeric(gmapRes[ , colId[i]]) >= minCoverage) & (as.numeric(gmapRes[ , (colId[i] + 1)]) >= minIdentity))
  
  curFilteredInfo <- matrix("NA", nrow = length(idx), 2)
  curFilteredInfo[ , 1] <- rownames(gmapRes)[idx]
  curFilteredInfo[ , 2] <- rep(colnames(gmapRes)[colId[i]], length(idx))
  
  filteredTransInfo <- as.matrix(rbind(filteredTransInfo, curFilteredInfo))
  
}

filteredTransInfo <- filteredTransInfo[-1, ]
colnames(filteredTransInfo) <- c("transcript_ID", "refSeq_ID")
filteredTransID <- unique(filteredTransInfo[ , 1])

write.table(filteredTransInfo, options$output1, sep = "\t", quote = F, col.names = T, row.names = F)
write.table(filteredTransID, options$output2, sep = "\t", quote = F, col.names = F, row.names = F)

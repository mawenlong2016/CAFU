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
option_specification = matrix(c('input', 'i', 2, 'character',
                                'output', 'o', 2, 'character'
                               ),
                              byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

FPKM <- as.matrix(read.table(file = options$input, sep = "", quote = "", header = T))
class(FPKM) <- "numeric"
idx <- apply(FPKM, 1, FUN = function(x) length(which(x>=1)))
idx <- which(idx > 0.1 * ncol(FPKM))

expTransID <- rownames(FPKM)[idx]

write.table(expTransID, file = options$output, sep = "\t", quote = F, col.names = F, row.names = F)

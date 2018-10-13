#### Command to run tool:
# Rscript /home/malab10/software/galaxy/tools/myTools/length_distribution.R --input1 $expMat --input2 $maxCutOff --input3 $minCutOff --output1 $transGC --outpout2 $transGCPlot;

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
                                'output1', 'o1', 2, 'character',
                                'output2', 'o2', 2, 'character'
                               ),
                              byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

library(ineq)
library(ggplot2)

allFPKM <- as.matrix(read.table(file = options$input1, sep = "\t", quote = "", header = T, row.names = 1))
allFPKM <- allFPKM[grep("Contig", rownames(allFPKM)), ]

class(allFPKM) <- "numeric"
allFPKM <- log2(allFPKM + 1)
idx <- which(apply(allFPKM, 1, sum) == 0)
if(length(idx) != 0){
  allFPKM <- allFPKM[-idx, ]
}

GC <- as.matrix(apply(allFPKM, 1, ineq))

write.table(GC, file = options$output1, sep = "\t", quote = F, col.names = F, row.names = T)

a <- c(1:nrow(GC))
b <- as.numeric(GC[ , 1])

GCdf <- data.frame(a, b)
pdf(file = options$output2, width = 5, height = 5)
ggplot(GCdf,
       aes(x = a,
           y = b)
       ) + 
  geom_point() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(name = "sample number"
                     ) +
  scale_y_continuous(name = "Gini coefficient"
                     )
dev.off()

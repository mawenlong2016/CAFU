#### Command to run tool:
#Rscript /home/chenzhuod/galaxy/tools/CAFU/10_length_distribution.R --input1 $novel_len --input2 $ref_len --output $plot_transcripts_length ;

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
                                'output', 'o', 2, 'character'
                               ),
                              byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

library(ggplot2)

noveLen <- as.matrix(read.table(file = options$input1, sep = "\t", quote = "", header = F))
refLen <- as.matrix(read.table(file = options$input2, sep = "\t", quote = "", header = F))

lenStat <- matrix("NA", nrow = 8, 2)
rownames(lenStat) <- c("1-300", "301-500", "501-1000", "1001-1500", "1501-2000", "2001-3000", "3001-4000", "4001-")
colnames(lenStat) <- c("novelTrans", "refTrans")

for( i in 1:length(rownames(lenStat))){
  if( i <= 7 ){
    novlIdx <- which((as.numeric(noveLen[ , 2]) >= as.numeric(strsplit(rownames(lenStat)[i], "-")[[1]][1])) & (as.numeric(noveLen[ , 2]) <= as.numeric(strsplit(rownames(lenStat)[i], "-")[[1]][2])))
    refIdx <- which((as.numeric(refLen[ , 2]) >= as.numeric(strsplit(rownames(lenStat)[i], "-")[[1]][1])) & (as.numeric(refLen[ , 2]) <= as.numeric(strsplit(rownames(lenStat)[i], "-")[[1]][2])))
    lenStat[i, 1] <- length(novlIdx)
    lenStat[i, 2] <- length(refIdx)
  } else {
    novlIdx <- which(as.numeric(noveLen[ , 2]) >= as.numeric(strsplit(rownames(lenStat)[i], "-")[[1]][1]))
    refIdx <- which(as.numeric(refLen[ , 2]) >= as.numeric(strsplit(rownames(lenStat)[i], "-")[[1]][1]))
    lenStat[i, 1] <- length(novlIdx)
    lenStat[i, 2] <- length(refIdx)
  }
}

lenStat[ , 1] <- (as.numeric(lenStat[ , 1]) / sum(as.numeric(lenStat[ , 1]))) * 100
lenStat[ , 2] <- (as.numeric(lenStat[ , 2]) / sum(as.numeric(lenStat[ , 2]))) * 100

pVal <- ks.test(as.numeric(lenStat[ , 1]), as.numeric(lenStat[ , 2]))$p.value

a <- rep(c(1:8), 2)
b <- as.numeric(lenStat)
c <- rep(colnames(lenStat), each = length(rownames(lenStat)))

lenStatDF <- data.frame(a, b, c)

pdf(file = options$output, width = 7, height = 5)
ggplot(lenStatDF,
       aes(x = a,
           y = b,
           fill = c)
) +
  geom_bar(stat = "identity",
           position = position_dodge(0.7),
           width = 0.7
  ) +
  scale_fill_manual(values = c("#1E499D", "#e71219")
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(name = paste0("length distribution (p-value = ", pVal, ")"),
                     breaks = 1:8,
                     labels = c("1-300", "301-500", "501-1000", "1001-1500", "1501-2000", "2001-3000", "3001-4000", "4001-")
  ) +
  theme(axis.text.x = element_text(vjust = 0.5,
                                   angle = 45)
  ) + 
  scale_y_continuous(name = "percent (%)")
dev.off()

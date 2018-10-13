#### Command to run tool:
# Rscript /home/malab10/software/galaxy/tools/myTools/GC_distribution.R --input1 $output_knownTran --input2 $output_novelTran --output $plot_GC_content

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

knownGC <- as.matrix(read.table(file = options$input1, sep = "\t", quote = "", header = F))
novelGC <- as.matrix(read.table(file = options$input2, sep = "\t", quote = "", header = F))

GCStat <- matrix("NA", 20, 2)
rownames(GCStat) <- c("0-5",  "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80", "80-85", "85-90", "90-95", "95-100")
colnames(GCStat) <- c("nvoelTransGC", "knownTransGC")

for( i in 1:length(rownames(GCStat))){
  
  novelIdx <- which((as.numeric(novelGC[ , 1]) >= (as.numeric(strsplit(rownames(GCStat)[i], "-")[[1]][1]) * 0.01)) & (as.numeric(novelGC[ , 1]) < (as.numeric(strsplit(rownames(GCStat)[i], "-")[[1]][2]) * 0.01)))
  knownIdx <- which((as.numeric(knownGC[ , 1]) >= (as.numeric(strsplit(rownames(GCStat)[i], "-")[[1]][1]) * 0.01)) & (as.numeric(knownGC[ , 1]) < (as.numeric(strsplit(rownames(GCStat)[i], "-")[[1]][2]) * 0.01)))
  
  GCStat[i, 1] <- length(novelIdx)
  GCStat[i, 2] <- length(knownIdx)
  
}

GCStat[ , 1] <- (as.numeric(GCStat[ , 1]) / sum(as.numeric(GCStat[ , 1]))) * 100
GCStat[ , 2] <- (as.numeric(GCStat[ , 2]) / sum(as.numeric(GCStat[ , 2]))) * 100

pVal <- ks.test(as.numeric(GCStat[ , 1]), as.numeric(GCStat[ , 2]))$p.value

a <- rep(c(1:20), 2)
b <- as.numeric(GCStat)
c <- rep(colnames(GCStat), each = length(rownames(GCStat)))

GCStatDF <- data.frame(a, b, c)

pdf(file = options$output, width = 7, height = 5)
ggplot(GCStatDF,
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
  scale_x_continuous(name = paste0("GC content distribution (p-value = ", pVal, ")"),
                     breaks = 1:20,
                     labels = c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80", "80-85", "85-90", "90-95", "95-100")
  ) +
  theme(axis.text.x = element_text(vjust = 0.5,
                                   angle = 45)
        ) + 
  scale_y_continuous(name = "percent (%)")
dev.off()

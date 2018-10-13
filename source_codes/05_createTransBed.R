setwd("/home/chenzhuod/galaxy/database/files/Expression/expression_evidence/")

transLen <- as.matrix(read.table("assembledTrans_len", sep = "\t", quote	= "", header = F))
for( i in 1:nrow(transLen) ){
	
	curBed <- matrix("0", 1, 3)
	curBed[1, 1] <- transLen[i, 1]
	curBed[1, 3] <- as.numeric(transLen[i, 2])
	
	outDir <- paste0("./transBed/", transLen[i, 1], ".bed")
	write.table(curBed, outDir, sep = "\t", quote = F, col.names = F, row.names = F)

}

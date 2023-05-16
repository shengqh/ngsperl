require(tidyverse)
dnmtools_path <- paste("..", "..", "dnmtools", "result", sep="/")

#prepare to find all specific format files
cpg.pat <- ".meth$"
cpg.all <- list.files(path = dnmtools_path, pattern = cpg.pat, all.files = FALSE, recursive = T, ignore.case = FALSE, include.dirs = FALSE)
cpg.all <- cpg.all[!(str_detect(cpg.all,"all.meth$"))]

#count average bvalue for each hmr of each sample
for (i in 1:length(cpg.all)) {
	#read in *.meth files
	cpg.infile <- paste(dnmtools_path, cpg.all[i], sep = "/")
	id <- gsub("(.*).meth", "\\1", cpg.all[i])
	#format cpg data frame
	DT.cpg <- read.table(cpg.infile, header = F, stringsAsFactors = FALSE)
	colnames(DT.cpg) <- c("chr", "base", "strand", "Type", "freqC", "coverage")
	DT.cpg$chrBase <- paste(DT.cpg$chr, DT.cpg$base, sep = ".")
	DT.cpg$freqC <- round(DT.cpg$freqC * 100, 3)
	DT.cpg$freqT <- 100 - DT.cpg$freqC
	DT.cpg$strand <- "F"
	DT.cpg <- DT.cpg[,c("chrBase","chr","base","strand","coverage","freqC","freqT")]
	write.table(DT.cpg, file = paste0("../result/", paste(id, "CpG", "txt", sep = ".")), sep = "\t", quote = F, row.names = F)
}


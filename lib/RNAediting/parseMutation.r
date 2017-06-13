options(bitmapType='cairo')

library(ggplot2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
inputfile = args[1]
outputprefix = args[2]
width=as.numeric(args[3])
height=as.numeric(args[4])
thresholds=as.numeric(args[5:length(args)])

#inputfile="H:/shengquanhu/projects/GuoYan/20160922_rnaediting/rnaediting_result/result/rnaediting.SNV.tsv"
#outputprefix="H:/shengquanhu/projects/GuoYan/20160922_rnaediting/rnaediting_result/result/rnaediting.SNV"
#width=6000
#height=3000
#thresholds=c(1, 2, 5)

snv=read.table(inputfile, header=T, sep="\t")

snv$RefAllele="REF"

for (threshold in thresholds) {
	cat("threshold = ", threshold, "\n")
	filtered = snv[snv$AltAllelePercentage >= threshold,]
	
	slim1=filtered[,c("File", "Position", "RefAllele", "RefAllelePercentage")]
	colnames(slim1)=c("File", "Position", "Allele", "AllelePercentage")
	slim2=filtered[,c("File", "Position", "AltAllele", "AltAllelePercentage")]
	colnames(slim2)=c("File", "Position", "Allele", "AllelePercentage")
	slim=rbind(slim1, slim2)
	slim$Allele=factor(slim$Allele, levels = c("REF", "A", "T", "G", "C"))
	
	myColors <- c("gray", "red", "blue", "green", "purple")
	names(myColors) <- levels(slim$Allele)
	
	outputfile = paste0(outputprefix, ".perc" , threshold, ".png")
	png(filename=outputfile, width=width, height=height, res=300)
	print(ggplot(data=slim, aes(x=Position, y=AllelePercentage)) + 
					geom_bar(stat = "identity", position="stack", aes(fill=Allele)) + 
					scale_fill_manual(values = myColors)+
					facet_grid(File~.) +
					theme_cowplot() +
					theme(strip.background=element_blank(),
							axis.ticks.y=element_blank(),
							axis.text.y=element_blank(),
							axis.line.y=element_blank(),
							strip.text.y = element_text(angle=0,hjust = 0),
							panel.background = element_rect(fill="gray95")
					) +
					ylab("Allele Percentage"))
	dev.off()
}


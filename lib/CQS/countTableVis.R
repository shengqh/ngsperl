#############################
#Vis for all count tables: Group 1, 2, 4; small Rna Category; rRNA categry; tRNA category
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
mappingResultFile<-parFile1
databaseLogFile<-parFile2
totalCountFile<-parFile3
groupFileList<-parSampleFile1

mappingResult<-read.delim(mappingResultFile,header=T,row.names=1, check.names=F)

if (databaseLogFile!="") { #Group Count Table
	databaseLog<-read.delim(databaseLogFile,header=T,as.is=T)
	id2Species<-databaseLog$Species
	names(id2Species)<-databaseLog$Id
	
	mappingResultExpand<-expandCountTableByName(mappingResult)
	speciesInMappingResult<-id2Species[row.names(mappingResultExpand)]
	
	mappingResult2Species<-aggregateCountTable(mappingResultExpand,speciesInMappingResult)

	#short name
	row.names(mappingResult2Species)<-sapply(strsplit(row.names(mappingResult2Species),"_"),function(x) {
				if (length(x)<=3) {
					paste(x,collapse="_")
				} else if (grepl("^\\d+$",x[2])) {
					paste(x[1:3],collapse="_")
				} else {
					paste(x[1:2],collapse="_")
				}
			})
	write.csv(mappingResult2Species,paste0(resultFile,".Species.csv"))
} else {
	mappingResult2Species<-mappingResult
}

#Pie chart for all samples
png(paste0(resultFile,".Piechart.png"),width=3000,height=3000,res=300)
p<-ggpie(mappingResult2Species,fill="Category", y="Reads",facet="Sample",maxCategory=maxCategory,textSize=textSize)
print(p)
dev.off()

#Barplot for all samples
if (totalCountFile!="") { #normlize with total count *10^6
	totalCount<-read.csv(totalCountFile,header=T,as.is=T,row.names=1)
	totalCount<-unlist(totalCount["Reads for Mapping",])
	mappingResult2Species<-10^6*t(t(mappingResult2Species)/totalCount[colnames(mappingResult2Species)])
}
width<-max(3000,75*ncol(mappingResult2Species))
#height<-max(1500,50*min(maxCategory+ncol(mappingResult2Species),nrow(mappingResult2Species)))
height<-1500
png(paste0(resultFile,".Barplot.png"),width=width,height=height,res=300)
tableBarplot(mappingResult2Species,maxCategory=maxCategory,textSize=textSize,ylab="Mapped Reads per Million")
dev.off()

#Pie chart for all groups
if (groupFileList!="") {
	sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)
	#keep the groups with samples in the count table
	sampleToGroup<-sampleToGroup[which(sampleToGroup[,1] %in% colnames(mappingResult)),]
	
	mappingResult2SpeciesBySampleGroup<-mergeTableBySampleGroup(mappingResult2Species,sampleToGroup)
	
	write.csv(mappingResult2SpeciesBySampleGroup,paste0(resultFile,".PercentGroups.csv"))
	
	png(paste0(resultFile,".Piechart.Group.png"),width=3000,height=3000,res=300)
	p<-ggpie(mappingResult2SpeciesBySampleGroup,fill="Category", y="Reads",facet="Sample",maxCategory=maxCategory,textSize=textSize)
	print(p)
	dev.off()
}


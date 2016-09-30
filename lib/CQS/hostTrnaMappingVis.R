#############################
#Vis for tRNA category
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
mappingResultFile<-parFile1
#databaseLogFile<-parFile2
totalCountFile<-parFile3
groupFileList<-parSampleFile1
groupVisLayoutFileList<-parSampleFile2

#tRNA count file and Tables
mappingResult<-read.delim(mappingResultFile,header=T,row.names=1, check.names=F)
isDataNumeric = unlist(lapply(mappingResult[1,], function(x){is.numeric(x)}))
if (any(isDataNumeric)) {
	index = 1
	while(!all(isDataNumeric[index:ncol(mappingResult)])){
		index = index + 1
	}
} else {
	cat("Error: No numeric data found for tRNA category figure \n")
	quit(save="yes")
}
mappingResult<-mappingResult[,c(index:ncol(mappingResult))]

mappingResultExpand<-expandCountTableByName(mappingResult)
#Make one table into different sub Tables, need to match tRNA name
if (grepl("^tRNA-[A-Za-z]+-[A-Za-z]+",row.names(mappingResultExpand)[1])) { #tRNA like tRNA-Gly-GCC-1-2
	nameSubtRNA<-sapply(strsplit(row.names(mappingResultExpand),"-"),function(x) paste(x[2:3],collapse=""))
	nameSubtRNA1<-sapply(strsplit(row.names(mappingResultExpand),"-"),function(x) x[2])
} else { #tRNA like chr6.tRNA1014-CysGCA
	nameLength<-nchar(row.names(mappingResultExpand))
	nameSubtRNA<-substr(row.names(mappingResultExpand),nameLength-5,nameLength)
	nameSubtRNA1<-substr(row.names(mappingResultExpand),nameLength-5,nameLength-3)
}
trnaCountTableExpandByRNA<-aggregateCountTable(mappingResultExpand,nameSubtRNA)
trnaCountTableExpandByRNA1<-aggregateCountTable(mappingResultExpand,nameSubtRNA1)

#Pie chart for tables
ggpieToFile(trnaCountTableExpandByRNA,fileName=paste0(resultFile,".tRNAType2.Piechart.png"),maxCategory=maxCategory,textSize=textSize)
ggpieToFile(trnaCountTableExpandByRNA1,fileName=paste0(resultFile,".tRNAType1.Piechart.png"),maxCategory=maxCategory,textSize=textSize)

#Barplot for tables
tableBarplotToFile(dat=trnaCountTableExpandByRNA,fileName=paste0(resultFile,".tRNAType2.Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)
tableBarplotToFile(dat=trnaCountTableExpandByRNA1,fileName=paste0(resultFile,".tRNAType1.Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart for tables
ggpieGroupToFile(dat=trnaCountTableExpandByRNA,fileName=paste0(resultFile,".tRNAType2.Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".tRNAType2.PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)
ggpieGroupToFile(dat=trnaCountTableExpandByRNA1,fileName=paste0(resultFile,".tRNAType1.Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".tRNAType1.PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)




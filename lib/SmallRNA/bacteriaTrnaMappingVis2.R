
source("countTableVisFunctions.R")

options(bitmapType='cairo')

#############################
#Vis for tRNA category
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
mappingResultFile<-parFile1
databaseLogFile<-parFile2
totalCountFile<-parFile3
groupFileList<-parSampleFile1
groupVisLayoutFileList<-parSampleFile2

facetColCount=getFacetColCount(groupFileList)

#Category count Table
mappingResultCategory<-read.delim(gsub(".count$",".Species.count",mappingResultFile),header=T,row.names=1, check.names=F)

#tRNA count file and Tables
mappingResult<-read.delim(mappingResultFile,header=T,row.names=1, check.names=F)
mappingResultExpand<-expandCountTableByName(mappingResult)

#Make one table into different sub Tables
nameSub<-strsplit(row.names(mappingResultExpand),"_tRNA-|.trna\\d+-")
temp<-sapply(nameSub,function(x) x[2])
temp<-gsub("-\\d+-\\d+$","",temp)
nameSubtRNA<-paste0(substr(temp,0,3),"-",substr(temp,nchar(temp)-2,nchar(temp)),sep="")
nameSubtRNA1<-substr(temp,0,3)

cats<-read.table(parFile2, sep="\t", header=T, comment='$')
if ("Species" %in% colnames(cats)){
	cats$Genome12<-sapply(strsplit(cats$Id,"_"),function(x) paste0(x[1:2],collapse="_"))
	cat1Map<-split(cats$Genome12, cats$Id)
	cats2<-cats[!duplicated(cats$Genome12),]
	cat2Map<-split(cats2$Species, cats2$Genome12)
	nameSubSpecies12<-unlist(cat1Map[row.names(mappingResultExpand)])
}else{
	cats$Genome12<-sapply(strsplit(cats$Genome," "),function(x) paste0(x[1:2],collapse="_"))
	cat1Map<-split(cats$Genome12, cats$Chrom)
	cats2<-cats[!duplicated(cats$Genome12),]
	cat2Map<-split(cats2$Domain, cats2$Genome12)
	nameSubSpecies12<-unlist(cat1Map[row.names(mappingResultExpand)])
}

#trnaCountTableExpandBySpecies<-aggregateCountTable(mappingResultExpand,nameSubSpecies)
trnaCountTableExpandBySpecies12<-aggregateCountTable(mappingResultExpand,nameSubSpecies12)
trnaCountTableExpandByRNA<-aggregateCountTable(mappingResultExpand,nameSubtRNA)
trnaCountTableExpandByRNA1<-aggregateCountTable(mappingResultExpand,nameSubtRNA1)

species12<-cbind(data.frame("Domain" = unlist(cat2Map[rownames(trnaCountTableExpandBySpecies12)])), trnaCountTableExpandBySpecies12)

write.csv(species12,paste0(resultFile,".Species12.csv"))
write.csv(trnaCountTableExpandByRNA,paste0(resultFile,".tRNAType2.csv"))
write.csv(trnaCountTableExpandByRNA1,paste0(resultFile,".tRNAType1.csv"))

#Make Individual Graphics
for ( i in 1:ncol(mappingResult)) {
	png(paste0(resultFile,".",colnames(mappingResult)[i],".tRNA.png"),height=3000,width=3000,res=300)
	par(mfrow=c(2,2))
	par(mar=c(2,2,2,2))
	
	temp<-as.matrix(trnaCountTableExpandBySpecies12)
	basicPie(temp[,i])
	
	temp<-as.matrix(trnaCountTableExpandByRNA)
	basicPie(temp[,i])
	
	temp<-as.matrix(trnaCountTableExpandByRNA1)
	basicPie(temp[,i])
	
	temp<-as.matrix(mappingResultCategory)
	basicPie(temp[,i])
	
	dev.off()
}


maxCategoryForSpecies12<-min(c(maxCategory,3),na.rm=T)

bacteriaSpecies12 = trnaCountTableExpandBySpecies12[species12$Domain %in% c("Bacteria","bacteria"),]
bacteriaResultFile = paste0(resultFile, ".Bacteria")

#Pie chart for tables
ggpieToFile(trnaCountTableExpandBySpecies12,fileName=paste0(resultFile,".Species12.Piechart.png"),maxCategory=maxCategoryForSpecies12,textSize=textSize,facetColCount=facetColCount)
ggpieToFile(bacteriaSpecies12,fileName=paste0(bacteriaResultFile,".Species12.Piechart.png"),maxCategory=maxCategoryForSpecies12,textSize=textSize,facetColCount=facetColCount)

ggpieToFile(trnaCountTableExpandByRNA,fileName=paste0(resultFile,".tRNAType2.Piechart.png"),maxCategory=maxCategory,textSize=textSize,facetColCount=facetColCount)
ggpieToFile(trnaCountTableExpandByRNA1,fileName=paste0(resultFile,".tRNAType1.Piechart.png"),maxCategory=maxCategory,textSize=textSize,facetColCount=facetColCount)
ggpieToFile(mappingResultCategory,fileName=paste0(resultFile,".Category.Piechart.png"),maxCategory=NA,textSize=textSize,facetColCount=facetColCount)

#Barplot for tables
tableBarplotToFile(dat=trnaCountTableExpandBySpecies12,fileName=paste0(resultFile,".Species12.Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategoryForSpecies12,textSize=textSize)
tableBarplotToFile(dat=bacteriaSpecies12,fileName=paste0(bacteriaResultFile,".Species12.Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategoryForSpecies12,textSize=textSize)

tableBarplotToFile(dat=trnaCountTableExpandByRNA,fileName=paste0(resultFile,".tRNAType2.Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)
tableBarplotToFile(dat=trnaCountTableExpandByRNA1,fileName=paste0(resultFile,".tRNAType1.Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)
tableBarplotToFile(dat=mappingResultCategory,fileName=paste0(resultFile,".Category.Barplot.png"),totalCountFile=totalCountFile,maxCategory=NA,textSize=textSize)

#Barplot for Group samples
tableBarplotToFile(dat=trnaCountTableExpandBySpecies12,fileName=paste0(resultFile,".Species12.Group.Barplot.png"),groupFileList=groupFileList,outFileName=paste0(resultFile,".Species12.ReadsPerMillionGroups.csv"),totalCountFile=totalCountFile,maxCategory=maxCategoryForSpecies12,textSize=textSize)
tableBarplotToFile(dat=bacteriaSpecies12,fileName=paste0(bacteriaResultFile,".Species12.Group.Barplot.png"),groupFileList=groupFileList,outFileName=paste0(resultFile,".Species12.ReadsPerMillionGroups.csv"),totalCountFile=totalCountFile,maxCategory=maxCategoryForSpecies12,textSize=textSize)

tableBarplotToFile(dat=trnaCountTableExpandByRNA,fileName=paste0(resultFile,".tRNAType2.Group.Barplot.png"),groupFileList=groupFileList,outFileName=paste0(resultFile,".tRNAType2.ReadsPerMillionGroups.csv"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)
tableBarplotToFile(dat=trnaCountTableExpandByRNA1,fileName=paste0(resultFile,".tRNAType1.Group.Barplot.png"),groupFileList=groupFileList,outFileName=paste0(resultFile,".tRNAType1.ReadsPerMillionGroups.csv"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)
tableBarplotToFile(dat=mappingResultCategory,fileName=paste0(resultFile,".Category.Group.Barplot.png"),groupFileList=groupFileList,outFileName=paste0(resultFile,".Category.ReadsPerMillionGroups.csv"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart for tables
ggpieGroupToFile(dat=trnaCountTableExpandBySpecies12,fileName=paste0(resultFile,".Species12.Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".Species12.PercentGroups.csv"),maxCategory=maxCategoryForSpecies12,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)
ggpieGroupToFile(dat=bacteriaSpecies12,fileName=paste0(bacteriaResultFile,".Species12.Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(bacteriaResultFile,".Species12.PercentGroups.csv"),maxCategory=maxCategoryForSpecies12,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)

ggpieGroupToFile(dat=trnaCountTableExpandByRNA,fileName=paste0(resultFile,".tRNAType2.Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".tRNAType2.PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)
ggpieGroupToFile(dat=trnaCountTableExpandByRNA1,fileName=paste0(resultFile,".tRNAType1.Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".tRNAType1.PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)
ggpieGroupToFile(dat=mappingResultCategory,fileName=paste0(resultFile,".Category.Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".Category.PercentGroups.csv"),maxCategory=NA,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)


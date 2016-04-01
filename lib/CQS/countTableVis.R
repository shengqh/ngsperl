#############################
#Vis for all count tables: Group 1, 2, 4; rRNA categry;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
mappingResultFile<-parFile1
databaseLogFile<-parFile2
totalCountFile<-parFile3
groupFileList<-parSampleFile1
groupVisLayoutFileList<-parSampleFile2

mappingResult<-read.delim(mappingResultFile,header=T,row.names=1, check.names=F)
mappingResult2Species<-countTableToSpecies(dat=mappingResult,databaseLogFile=databaseLogFile,outFileName=paste0(resultFile,".Species.csv"))

#Pie chart for all samples
ggpieToFile(mappingResult2Species,fileName=paste0(resultFile,".Piechart.png"),maxCategory=maxCategory,textSize=textSize)

#Barplot for all samples
tableBarplotToFile(dat=mappingResult2Species,fileName=paste0(resultFile,".Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart
ggpieGroupToFile(dat=mappingResult2Species,fileName=paste0(resultFile,".Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".PercentGroups.csv"),maxCategory=maxCategory,textSize=textSize,visLayoutFileList=groupVisLayoutFileList)

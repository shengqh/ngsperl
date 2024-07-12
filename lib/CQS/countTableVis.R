rm(list=ls()) 
outFile='DM_7369.fungus_group4Mapping.Result'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parFile1='/nobackup/vickers_lab/projects/20240613_7369_DM_smRNA_human_byMars_rerun/nonhost_genome/bowtie1_fungus_group4_pm_table/result/fungus_group4_pm_DM_7369.Species.count'
parFile2=''
parFile3='/nobackup/vickers_lab/projects/20240613_7369_DM_smRNA_human_byMars_rerun/preprocessing/fastqc_post_trim_summary/result/DM_7369.countInFastQcVis.Result.Reads.csv'
maxCategory=4;

setwd('/nobackup/vickers_lab/projects/20240613_7369_DM_smRNA_human_byMars_rerun/data_visualization/nonhost_genome_fungus_group4_vis/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
options(bitmapType='cairo')

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

myoptions=read_file_map(parSampleFile3, do_unlist=FALSE)
barplot_width_px=as.integer(myoptions$barplot_width_px)
barplot_height_px=as.integer(myoptions$barplot_height_px)
textSize=as.integer(myoptions$textSize)
groupTextSize=as.integer(myoptions$groupTextSize)

facetColCount=getFacetColCount(groupFileList)

mappingResult<-read.delim(mappingResultFile,header=T,row.names=1, check.names=F)
mappingResult2Species<-countTableToSpecies(dat=mappingResult,databaseLogFile=databaseLogFile,outFileName=paste0(resultFile,".Species.csv"))

stringLengthCut=40
if (max(nchar(row.names(mappingResult2Species)))>stringLengthCut) {
	row.names(mappingResult2Species)=make.unique(stringr::str_trunc(row.names(mappingResult2Species), stringLengthCut))
}

#Pie chart for all samples
ggpieToFile(mappingResult2Species,
            fileName=paste0(resultFile,".Piechart.png"),
            maxCategory=maxCategory,
            textSize=textSize,
            facetColCount=facetColCount)

#Barplot for all samples
tableBarplotToFile( dat=mappingResult2Species, 
                    fileName=paste0(resultFile,".Barplot.png"),
                    totalCountFile=totalCountFile,
                    maxCategory=maxCategory,
                    textSize=textSize,
                    xlab="",
                    width=barplot_width_px, 
                    height=barplot_height_px)

#Barplot for Group samples
tableBarplotToFile( dat=mappingResult2Species,
                    fileName=paste0(resultFile,".Group.Barplot.png"),
                    groupFileList=groupFileList,
                    outFileName=paste0(resultFile,".ReadsPerMillionGroups.csv"),
                    totalCountFile=totalCountFile,
                    maxCategory=maxCategory,
                    textSize=textSize,
                    xlab="",
                    height=barplot_height_px)

#Group Pie chart
ggpieGroupToFile( dat=mappingResult2Species,
                  fileName=paste0(resultFile,".Group.Piechart.png"),
                  groupFileList=groupFileList,
		              outFileName=paste0(resultFile,".PercentGroups.csv"),
                  maxCategory=maxCategory,
                  textSize=groupTextSize,
                  visLayoutFileList=groupVisLayoutFileList)

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')

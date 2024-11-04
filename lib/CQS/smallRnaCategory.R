rm(list=ls()) 
outFile='nextflex'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1=''
parFile2=''
parFile3=''
textSize=12;groupTextSize=12;

setwd('/nobackup/vickers_lab/projects/20240201_smallRNA_nextflex_comparison_hg38_byTiger/host_genome/bowtie1_genome_1mm_NTA_smallRNA_category/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
options(bitmapType='cairo')

option_tbl=fread(parSampleFile4, header=FALSE, data.table=FALSE) |> dplyr::rename(option=2, value=1)
myoptions=split(option_tbl$value, option_tbl$option)
color_theme=myoptions$color_theme

host_reads_colors=get_host_reads_colors(color_theme)
host_smallRNA_colors=get_host_smallRNA_colors(color_theme)

taskName<-outFile
categoryFileList<-parSampleFile1
groupFileList<-parSampleFile2
groupVisLayoutFileList<-parSampleFile3

if(!exists("visLayoutAlphabet")){
  visLayoutAlphabet<-FALSE
}

facetColCount=getFacetColCount(groupFileList)

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")

################################
#Generate category count table
################################
categoryFiles<-read.delim(categoryFileList,as.is=T,header=F)

categoryAll<-NULL
categoryFigure<-NULL
for (i in 1:nrow(categoryFiles)) {
  categoryFile<-categoryFiles[i,1]
  sampleName<-categoryFiles[i,2]
  categoryOne<-read.delim(categoryFile,as.is=T,header=F,comment.char = "#")
  categoryOne$Sample<-sampleName
  categoryAll<-rbind(categoryAll,categoryOne)
  
  row.names(categoryOne)<-categoryOne[,1]
  
  smallRna<-categoryOne["FeatureReads",2]
  genomes<-categoryOne["GenomeReads",2]
  tooshorts<-categoryOne["TooShortReads",2]
  unmapped<-categoryOne["UnannotatedReads",2]
  
  categoryOneFigure<-data.frame(Category=c("Unmapped","Too short", "Genome", "Small RNA"),
                                Count=c(unmapped,tooshorts, genomes,smallRna),
                                Sample=sampleName)
  categoryFigure<-rbind(categoryFigure,categoryOneFigure)
}
categoryFigure$Sample<-factor(categoryFigure$Sample,levels=sort(unique(as.character(categoryFigure$Sample))))

categoryAllTable<-acast(categoryAll,V1~Sample,value.var="V2")
summaryInd<-which(row.names(categoryAllTable) %in% c("TotalReads","MappedReads","FeatureReads","GenomeReads","TooShortReads","UnannotatedReads"))
categoryAllTable1<-categoryAllTable[summaryInd,,drop=F]
categoryAllTable2<-categoryAllTable[-summaryInd,,drop=F]
write.csv(rbind(categoryAllTable1,categoryAllTable2),paste0(taskName,".Category.Table.csv"))

################################
#Generate graphics for count tables
################################

if(!exists("drawInvidividual") || drawInvidividual){
  #Make individual figures
  for (i in 1:ncol(categoryAllTable1)) {
    png(paste0(colnames(categoryAllTable1)[i],".Category.png"),res=300,width=3000,height=1500)
    par(mfrow=c(1,2))
    par(mar=c(9,6,3,6))
    barplot(categoryAllTable1[c("TotalReads","MappedReads","FeatureReads", "TooShortReads"),i],
      names.arg=c("Total Reads","Mapped Reads","Small RNA", "Too Short Reads"),space=0.5,las=2,
      mar=c(9,9,9,9), col=rainbow(4))
    basicPie(categoryAllTable2[,i],addPercent=T)
   dev.off()
  }
}

#Pie Chart for Tables
p1<-ggpieToFile(categoryFigure,
                fileName=paste0(taskName,".Category1.Piechart.png"),
                maxCategory=NA,
                textSize=textSize,
                y="Count",
                transformTable=FALSE,
                colorNames=host_reads_colors,
                reOrder=FALSE,
                facetColCount=facetColCount)

p2<-ggpieToFile(categoryAllTable2,
                fileName=paste0(taskName,".Category2.Piechart.png"),
                maxCategory=NA,
                textSize=textSize,
                colorNames=host_smallRNA_colors, 
                facetColCount=facetColCount)

#Barplot for Tables
width<-max(2000,60*ncol(categoryAllTable2) + 200)

tableBarplotToFile( dat=categoryFigure,
                    fileName=paste0(taskName,".Category1.Barplot.png"),
                    totalCountFile="",
                    maxCategory=NA,
                    textSize=textSize,
                    transformTable=F,
                    colorNames=host_reads_colors,
                    y="Count", 
                    height=1200, 
                    width=width)

tableBarplotToFile( dat=categoryAllTable2,
                    fileName=paste0(taskName,".Category2.Barplot.png"),
                    totalCountFile="",
                    maxCategory=NA,
                    textSize=textSize, 
                    colorNames=host_smallRNA_colors, 
                    height=1200, 
                    width=width)

#Group Pie Chart for Tables
ggpieGroupToFile( dat=categoryFigure,
                  fileName=paste0(taskName,".Category1.Group.Piechart.png"),
                  groupFileList=groupFileList,
                  outFileName=paste0(taskName,".Category1.PercentGroup.Table.csv"),
                  maxCategory=NA,
                  textSize=groupTextSize,
                  y="Count",
                  transformTable=FALSE,
                  colorNames=host_reads_colors,
                  visLayoutFileList=groupVisLayoutFileList,
                  visLayoutAlphabet=visLayoutAlphabet)

ggpieGroupToFile( dat=categoryAllTable2,
                  fileName=paste0(taskName,".Category2.Group.Piechart.png"),
                  groupFileList=groupFileList,
                  outFileName=paste0(taskName,".Category2.PercentGroup.Table.csv"),
                  maxCategory=NA,
                  textSize=groupTextSize,
                  colorNames=host_smallRNA_colors, 
                  visLayoutFileList=groupVisLayoutFileList,
                  visLayoutAlphabet=visLayoutAlphabet)


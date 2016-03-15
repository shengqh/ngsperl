categoryFileList<-parSampleFile1
groupFileList<-parSampleFile2
taskName<-parFile1


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
	Unmapped<-categoryOne["TotalReads",2]-categoryOne["MappedReads",2]
	smallRna<-sum(categoryOne[-c(1:3),2])
	otherMapped<-categoryOne["MappedReads",2]-smallRna
	categoryOneFigure<-data.frame(Category=c("Unmapped","Other Mapped","Small RNA"),Count=c(Unmapped,otherMapped,smallRna),Sample=sampleName)
	categoryFigure<-rbind(categoryFigure,categoryOneFigure)
}

categoryAllTable<-acast(categoryAll,V1~Sample,value.var="V2")
write.csv(categoryAllTable,paste0(taskName,".Category.Table.csv"))

summaryInd<-which(row.names(categoryAllTable) %in% c("TotalReads","MappedReads","FeatureReads"))
categoryAllTable1<-categoryAllTable[summaryInd,]
categoryAllTable2<-categoryAllTable[-summaryInd,]

for (x in 1:ncol(categoryAllTable1)) {
	png(paste0(colnames(categoryAllTable1)[i],".Category.png"),res=300,width=3000,height=1500)
	par(mfrow=c(1,2))
	barplot(categoryAllTable1[c("TotalReads","MappedReads","FeatureReads"),i],
			names.arg=c("Total Reads","Mapped Reads","Small RNA"),space=0.5,las=2,
			mar=c(10,10,10,10), col=rainbow(3))
	tablePie(categoryAllTable2[,i],addPercent=T)
	dev.off()
}

if (groupFileList!="") {
	sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)
	#keep the groups with samples in the count table
	sampleToGroup<-sampleToGroup[which(sampleToGroup[,1] %in% colnames(categoryAllTable)),]
	
	categoryAllTable1Group<-mergeTableBySampleGroup(categoryAllTable1,sampleToGroup)
	categoryAllTable2Group<-mergeTableBySampleGroup(categoryAllTable2,sampleToGroup)
	
	colors<-makeColors(nrow(categoryAllTable2Group), "Set1")
	png(paste0(taskName,".Category.Group.Piechart.png"),width=2000,height=2000,res=300)
	p<-ggpie(categoryAllTable2Group,fill="Species", y="Reads",facet="Sample")
	print(
			p+scale_fill_manual(values=colors)
	)
	dev.off()
	
#	groupNames<-colnames(categoryAllTable2Group)
#	for (i in 1:length(groupNames)) {
#		png(paste0(groupNames[i],".Group.Category.png"),width=2000,height=1500,res=300)
#		par(mar=c(2,9,2,9))
#		temp<-as.matrix(categoryAllTable2Group)
#		groupPie(temp[,i],main=paste0("Group: ",groupNames[i]),addPercent=T)
#		dev.off()
#	}
}

colors<-makeColors(3, "Set1")
width<-max(3000,80*length(unique(categoryFigure$Sample)))
height<-max(1500,40*length(unique(categoryFigure$Sample)))
png(paste0(paste0(taskName,".Category.Barplot.png")), width=width, height=height, res=300)
p<-tableBarplot(categoryFigure,x="Sample",y="Count",fill="Category",transformTable=F)
print(p+scale_fill_manual(values=colors))
dev.off()











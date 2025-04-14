options(bitmapType='cairo')

projectName<-parFile1
projectDir<-parFile2

#projectTemplateFile<-"/home/zhaos/source/r_cqs/vickers/codesToPipeline/ProjectReport.Rmd"
#taskTemplateFile<-"/home/zhaos/source/r_cqs/vickers/codesToPipeline/TaskReport.Rmd"
projectTemplateFile<-"ProjectReport.Rmd"
taskTemplateFile<-"TaskReport.Rmd"

#colInProjectTable<-c("TaskName","SampleName","fileSize","Result")
colInTaskTable<-c("StepName","TaskName","SampleName","FileList","FileSize","Result")

projectReportDir<-paste0(projectDir,"/sequencetask/result")
resultFileListFile<-paste0(projectName,"_st_expect_result.tsv.check.csv")

library(rmarkdown)
library(DT)
expandDataTableByName<-function(x,sep=";") {
	namesEach<-strsplit(row.names(x),sep)
	namesEachUnlist<-unlist(namesEach)
	namesEachLength<-sapply(namesEach,length)
	result<-x[rep(seq.int(1,nrow(x)), namesEachLength),]
	row.names(result)<-namesEachUnlist
	colToExpand<-which(apply(x,2,function(y) any(grep(sep,y))))
	if (length(colToExpand)>0) {
		for (i in colToExpand) {
			result[,i]<-unlist(strsplit(x[,i],sep))
		}
	}
	return(result)
}
aggregateCountTable<-function(x,group,method=sum) {
	result<-aggregate(x, list(factor(group,levels=unique(group))), method)
	row.names(result)<-result[,1]
	result<-result[,-1]
	return(result)
}


addTag<-function(x,tag='a',attribute='href',label="",url="",urlBefore="",urlAfter="") {
	if (label[1]=="") {
		label<-x
	}
	if (url[1]=="") {
		url<-x
	}
	return(paste0('<',tag,' ',attribute,'="',urlBefore,url,urlAfter,'">',label,'</',tag,'>'))
}


#read result table
setwd(projectReportDir)
resultFileList<-read.csv(resultFileListFile,header=T,as.is=T)
#Only display files in result folder
temp<-grep("\\/result\\/.*$",resultFileList$FileList)
if (length(temp>0)) {
	resultFileList<-resultFileList[temp,]
}

#Prepare project report
stepNames<-unique(resultFileList$StepName)
figureToDisply<-data.frame(Title="Task Status",File=c(paste0(projectReportDir,"/",projectName,c(paste0("_st_expect_result.tsv_", stepNames, ".png"), paste0("_st_expect_result.tsv_", stepNames, ".RelativeFileSize.png")))),stringsAsFactors=F)

resultFileList$TaskResultFolder<-gsub( "\\/result\\/.*$", "\\/result/", resultFileList$FileList)
hasProjectDir<-grep(projectDir,resultFileList$TaskResultFolder)
if (length(hasProjectDir)>0) {
  resultFileList$TaskResultFolder[hasProjectDir]<-gsub(projectDir,"",resultFileList$TaskResultFolder[hasProjectDir])
  resultFileList$TaskResultFolder[hasProjectDir]<-paste0("../../",resultFileList$TaskResultFolder[hasProjectDir])
}
projectResultUnique<-aggregateCountTable(resultFileList,resultFileList$TaskName,function(x) x[1])

projectResultUniqueTable<-projectResultUnique[,c("StepName","TaskName","Result")]
colnames(projectResultUniqueTable)<-c("Step","Task","Status")
projectResultUniqueTable$'Result Folder'<-addTag(projectResultUnique$TaskResultFolder,label="Open Folder")
projectResultUniqueTable$Report<-addTag(projectResultUnique$TaskResultFolder,label="Open Report",urlAfter="TaskReport.html")
statusColor<-rep("black",nrow(projectResultUniqueTable))
statusColor[which(projectResultUniqueTable$Status=="WARN")]<-"skyblue"
statusColor[which(projectResultUniqueTable$Status=="FAIL")]<-"red"
projectResultUniqueTable$Status<-addTag(x=projectResultUniqueTable$Status,tag='font',attribute='color',url=statusColor)

#Render project report
try(render(projectTemplateFile,output_dir="."))

#Prepare and render task report
for (TaskFolder in unique(resultFileList$TaskResultFolder)) {
	#	resultFileListTask<-resultFileList[which(resultFileList$TaskResultFolder==TaskFolder),]
	TaskFolderKey<-gsub("\\.\\.\\/","",TaskFolder)
	temp<-grep(TaskFolderKey,resultFileList$FileList)
	if (length(temp)>0) {
		resultFileListTask<-resultFileList[temp,]
	} else {
		next;
	}
	
	row.names(resultFileListTask)<-resultFileListTask$FileList
	TaskName<-paste(unique(resultFileListTask$TaskName),collapse=", ")
	
	resultFileListTaskEachFile<-expandDataTableByName(resultFileListTask,sep=",")
	temp<-grep(TaskFolderKey,row.names(resultFileListTaskEachFile))
	resultFileListTaskEachFile<-resultFileListTaskEachFile[temp,]
	resultFileListTaskEachFile$FileList<-gsub( "^.*\\/result\\/", "", resultFileListTaskEachFile$FileList)
	
	resultFileListTaskEachFileTable<-resultFileListTaskEachFile[,colInTaskTable]
	resultFileListTaskEachFileTable$FileList<-addTag(resultFileListTaskEachFileTable$FileList)
	
	try(render(taskTemplateFile,output_dir=TaskFolder))
	
}


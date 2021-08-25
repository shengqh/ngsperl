options(bitmapType='cairo')

# TODO: Add comment
# 
# Author: zhaos
###############################################################################

#top100ReadFile<-"/scratch/cqs/zhaos/vickers/20160909_smallRNA_3018-KCV-77_78_79_mouse/class_independent/identical_sequence_count_table/result/KCV-77_78_79_sequence.read.count"
#readFiles<-c("/scratch/cqs/zhaos/vickers/20160909_smallRNA_3018-KCV-77_78_79_mouse/nonhost_library/bowtie1_tRNA_pm_table/result/tRNA_pm_KCV-77_78_79.read.count",
#		"/scratch/cqs/zhaos/vickers/20160909_smallRNA_3018-KCV-77_78_79_mouse/nonhost_library/bowtie1_rRNA_pm_table/result/rRNA_pm_KCV-77_78_79.read.count")
#readFilesModule<-c("Non-host tRNA","Non-host rRNA")
resultFile<-outFile
top100ReadFile<-parFile1
readFileList<-parSampleFile1
groupFileList<-parSampleFile2
groupVisLayoutFileList<-parSampleFile3

top100Table<-read.delim(top100ReadFile,header=T,row.names=1,as.is=T,check.names=F)
readFiles<-read.delim(readFileList,header=F,as.is=T)[,1]
mappingTable<-NULL
for (readFile in readFiles) {
	readsTable<-read.delim(readFile,header=T,row.names=1,as.is=T,check.names=F)
	mappingResult<-rep("N",nrow(top100Table))
	names(mappingResult)<-row.names(top100Table)
	
	temp<-intersect(row.names(top100Table),row.names(readsTable))
	mappingResult[temp]<-"Y"
	mappingTable<-cbind(mappingTable,mappingResult)
}
colnames(mappingTable)<-readFilesModule
temp<-apply(mappingTable,1,function(x) length(which(x=="Y")))
mappingTable<-cbind(mappingTable,NumberOfModulesMapped=temp)

result<-cbind(mappingTable,top100Table)
write.csv(result,paste0(resultFile,".ReadsMapping.Summary.csv"))


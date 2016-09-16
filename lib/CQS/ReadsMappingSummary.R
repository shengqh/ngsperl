# TODO: Add comment
# 
# Author: zhaos
###############################################################################

#top100ReadFile<-"/scratch/cqs/zhaos/vickers/20160909_smallRNA_3018-KCV-77_78_79_mouse/class_independent/identical_sequence_count_table/result/KCV-77_78_79_sequence.read.count"
#readFiles<-c("/scratch/cqs/zhaos/vickers/20160909_smallRNA_3018-KCV-77_78_79_mouse/nonhost_library/bowtie1_tRNA_pm_table/result/tRNA_pm_KCV-77_78_79.read.count",
#		"/scratch/cqs/zhaos/vickers/20160909_smallRNA_3018-KCV-77_78_79_mouse/nonhost_library/bowtie1_rRNA_pm_table/result/rRNA_pm_KCV-77_78_79.read.count")
#readFilesModule<-c("Non-host tRNA","Non-host rRNA")

top100Table<-read.delim(top100ReadFile,header=T,row.names=1,as.is=T)
mappingTable<-NULL
for (readFile in readFiles) {
	readsTable<-read.delim(readFile,header=T,row.names=1,as.is=T)
	mappingResult<-rep("N",nrow(top100Table))
	names(mappingResult)<-row.names(top100Table)
	
	temp<-intersect(row.names(top100Table),row.names(readsTable))
	mappingResult[temp]<-"Y"
	mappingTable<-cbind(mappingTable,mappingResult)
}
colnames(mappingTable)<-readFilesModule

result<-cbind(mappingTable,top100Table)



#parSampleFile1='fileList1.txt'
#outFile="20200729_UMI_somaticCNV_simpleProcessing"

cnvSegFiles = parSampleFile1

HasAnnotationFile=FALSE
if (file.exists(parSampleFile2)) {
  cnvSegAnnotationFiles=parSampleFile2
  cnvSegAnnotationFilesTable=read.delim(cnvSegAnnotationFiles,header=F,as.is=T,check.names=F)
  if (nrow(cnvSegAnnotationFilesTable)>0 & file.exists(cnvSegAnnotationFilesTable[1,1])) {
    HasAnnotationFile=TRUE
  }
}

minProbNum=5
chrToKeep=c(1:22,"X","Y",paste0("chr",c(1:22,"X","Y")))

# setwd("/scratch/cqs/zhaos/Pierre/WES/20200729_UMI_somaticCNV_simpleProcessing/files_somaticCNV_02_call/result/cromwell_finalOutputs")
# files=list.files(pattern=".hg38.called.seg$")
# samples=gsub(".hg38.called.seg$","",files)
# cnvSegFiles=data.frame(V1=paste0("/scratch/cqs/zhaos/Pierre/WES/20200729_UMI_somaticCNV_simpleProcessing/files_somaticCNV_02_call/result/cromwell_finalOutputs/",files),V2=samples)
# #remove PBMC/Normal
# cnvSegFiles=cnvSegFiles[-grep("PBMC",cnvSegFiles[,2]),]
# write.table(cnvSegFiles,"/scratch/cqs/zhaos/Pierre/WES/20200729_UMI_somaticCNV_simpleProcessing/somaticCNV_03_summary/fileList1.txt",row.names=FALSE,quote=FALSE,sep="\t",col.names=FALSE)

cnvSegFilesTable=read.delim(cnvSegFiles,header=F,as.is=T,check.names=F)

#find matched af segment files
afDir="/out/"
afFileNameSuffix=".af.igv.seg"

segFile=cnvSegFilesTable[,1]
afFileName=gsub(".called.seg",afFileNameSuffix,basename(segFile))
afFile=paste0(dirname(segFile),afDir,afFileName)
#afFile[which(!file.exists((afFile)))]=""
cnvSegFilesTable$AF=afFile

#process all segment files and merge them together
segmentAll=NULL
for (i in 1:nrow(cnvSegFilesTable)) {
  segmentFile=cnvSegFilesTable[i,1]
  segmentAfFile=cnvSegFilesTable[i,3]
  sampleName=cnvSegFilesTable[i,2]

  print(sampleName)
  segmentTableOne=read.delim(segmentFile,as.is=TRUE,header=TRUE,comment="@",check.names=F)
  segmentAfTableOne=read.delim(segmentAfFile,as.is=TRUE,header=TRUE,check.names=F)

  if (HasAnnotationFile) {
    segmentTableAnnotationOne=read.delim(cnvSegAnnotationFilesTable[i,1],as.is=TRUE,header=TRUE,check.names=F)
  }

  if (
    identical(segmentAfTableOne[,"Chromosome"],segmentTableOne[,"CONTIG"]) &
    identical(segmentAfTableOne[,"Start"],segmentTableOne[,"START"]) &
    identical(segmentAfTableOne[,"End"],segmentTableOne[,"END"])
  ) {
    #merge seg and af file, should start from seg file as af file may not exist
    segmentAllOneSample=data.frame(segmentAfTableOne,segmentTableOne[,-c(1:3)])
    if (HasAnnotationFile) {
      segmentAllOneSample=data.frame(segmentAllOneSample,segmentTableAnnotationOne[,c("start_gene","end_gene","genes","start_exon","end_exon","ref_allele","alt_allele")])
    }
    segmentAll=rbind(segmentAll,segmentAllOneSample)
  } else { #segements in seg and af files didn't match
    stop(paste0("segements in seg and af files didn't match"))
  }
}

segmentAll$copy_number=2^(segmentAll$MEAN_LOG2_COPY_RATIO)*2
chrYInd=which(segmentAll$Chromosome=="chrY" | segmentAll$Chromosome=="Y")
if (length(chrYInd)>0) { #chrY, normal copy_number=1, not 2
  segmentAll$copy_number[chrYInd]=2^(segmentAll$MEAN_LOG2_COPY_RATIO)
}
segmentAll$major_CN=(1-segmentAll$Segment_Mean)*segmentAll$copy_number
segmentAll$minor_CN=segmentAll$Segment_Mean*segmentAll$copy_number
segmentAll$segment_length=segmentAll$End-segmentAll$Start

#change column name to make it clear
colnames(segmentAll)=gsub("Segment_Mean","Segment_Mean_AF",colnames(segmentAll))
colnames(segmentAll)=gsub("Num_Probes","Num_Probes_AF",colnames(segmentAll))


#filtering: called !="0" and NUM_POINTS_COPY_RATIO>=minProbNum
#segmentAllFilter=segmentAll[which(segmentAll$Num_Probes_AF>minProbNum & segmentAll$NUM_POINTS_COPY_RATIO>minProbNum & segmentAll$Chromosome %in% chrToKeep),]
segmentAllFilter=segmentAll[which(segmentAll$CALL!="0" & segmentAll$NUM_POINTS_COPY_RATIO>=minProbNum & segmentAll$Chromosome %in% chrToKeep),]

write.table(segmentAll,paste0(outFile,".allCNV.seg"),row.names=FALSE,quote=FALSE,sep="\t")
write.table(segmentAllFilter,paste0(outFile,".allCNV.filter.seg"),row.names=FALSE,quote=FALSE,sep="\t")




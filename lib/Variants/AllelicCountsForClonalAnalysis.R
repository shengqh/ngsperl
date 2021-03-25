addAllelicCountsForClonalAnalysis

library(data.table)
library(GenomicRanges)

#mafFile to define mutation sites for all sample in clonal analysis
mafFile="/scratch/cqs/zhaos/workSync/Pierre/20200730_UMI_MuTect2PON_simpleProcessing/UMIMutect2PON.filter.allSamples.maf"
mafTable=fread(mafFile)

cnvFile="/scratch/cqs/zhaos/Pierre/WES/20200729_UMI_somaticCNV_simpleProcessing/somaticCNV_03_summary/20200729_UMI_somaticCNV_simpleProcessing.allCNV.filter.seg"
cnvTable=fread(cnvFile)

AllelicCountsDir="/scratch/cqs/zhaos/Pierre/WES/20200730_UMI_MuTect2PON_simpleProcessing/files_CollectAllelicCounts/result/cromwell_finalOutputs/"
AllelicCountsFiles=list.files(AllelicCountsDir,pattern=".allelicCounts.tsv")
AllelicCountsSampleToFile=data.frame(V1=paste0(AllelicCountsDir,AllelicCountsFiles),V2=gsub(".hg38.allelicCounts.tsv","",AllelicCountsFiles))


#patient level data
#sampleNames=c("Brushing_10658","Tumor_10658_A","Tumor_10658_B","Tumor_10658_C")
sampleNames=c("Tumor_10658_A","Tumor_10658_B","Tumor_10658_C")


mafTablePatientOne=mafTable[Tumor_Sample_Barcode %in% sampleNames & Variant_Type=="SNP",]
#posCol=c("Chromosome","Start_Position","End_Position")
#mafTablePatientOne[,..posCol]
#mafTablePatientOneGRange=makeGRangesFromDataFrame(mafTablePatientOne[,..posCol],seqnames.field=posCol[1],start.field=posCol[2],end.field=posCol[3],ignore.strand=TRUE)
#mafTablePatientOneGRange=reduce(mafTablePatientOneGRange)
posCol=c("Chromosome","Start_Position")
mafTablePatientOneRangeKey=mafTablePatientOne[,..posCol]
colnames(mafTablePatientOneRangeKey)=c("CONTIG","POSITION")
mafTablePatientOneRangeKey=unique(mafTablePatientOneRangeKey)
mafTablePatientOneRangeKey=mafTablePatientOneRangeKey[order(mafTablePatientOneRangeKey),]

#cnvTablePatientOne=cnvTable[Sample %in% sampleNames,]
#cnvTablePatientOneGRange=makeGRangesFromDataFrame(cnvTablePatientOne,ignore.strand=TRUE)


AllelicCountsTableListSciClone=list()
tableAllPyClone=NULL
cnvTableListSciClone=list()
for (sampleOne in sampleNames) {
  AllelicCountsFileOne=AllelicCountsSampleToFile[which(AllelicCountsSampleToFile[,2]==sampleOne),1]
  AllelicCountsTableOne=fread(AllelicCountsFileOne,skip="CONTIG")
  setkey(AllelicCountsTableOne,CONTIG,POSITION)
  AllelicCountsTableOne=AllelicCountsTableOne[mafTablePatientOneRangeKey,]
  
  AllelicCountsTableOneSciClone=data.frame(AllelicCountsTableOne[,c("CONTIG","POSITION","REF_COUNT","ALT_COUNT")])
  AllelicCountsTableOneSciClone$VAF=AllelicCountsTableOneSciClone$ALT_COUNT*100/(AllelicCountsTableOneSciClone$REF_COUNT+AllelicCountsTableOneSciClone$ALT_COUNT)
  AllelicCountsTableListSciClone=c(AllelicCountsTableListSciClone,list(AllelicCountsTableOneSciClone))
  
  cnvTableOneSciClone=data.frame(cnvTable[Sample==sampleOne,c("Chromosome","Start","End","copy_number")])
  cnvTableListSciClone=c(cnvTableListSciClone,list(cnvTableOneSciClone))
  
  ####################
  #PyClone-VI
  ####################
  AllelicCountsTableOnePyClone=AllelicCountsTableOne
  AllelicCountsTableOnePyClone$End=AllelicCountsTableOnePyClone$POSITION
  
  cnvTableOnePyClone=cnvTable[Sample==sampleOne,c("Chromosome","Start","End","major_CN","minor_CN")]
  setkey(cnvTableOnePyClone, Chromosome, Start, End)
  
  overlapResult=foverlaps(AllelicCountsTableOnePyClone, cnvTableOnePyClone, by.x=c("CONTIG", "POSITION", "End"),type="any", which=TRUE)
  
  tableOnePyClone=cbind(AllelicCountsTableOnePyClone,cnvTableOnePyClone[overlapResult$yid,])
  tableOnePyClone$mutation_id=paste0(tableOnePyClone$CONTIG,"_",tableOnePyClone$POSITION,"_",tableOnePyClone$REF_NUCLEOTIDE,"_",tableOnePyClone$ALT_NUCLEOTIDE)
  tableOnePyClone$sample_id=sampleOne
  tableOnePyClone$normal_cn=2
  
  tableOnePyClone=tableOnePyClone[,c("mutation_id","sample_id","REF_COUNT","ALT_COUNT","major_CN","minor_CN","normal_cn")]
  colnames(tableOnePyClone)=gsub("REF_COUNT","ref_counts",colnames(tableOnePyClone))
  colnames(tableOnePyClone)=gsub("ALT_COUNT","alt_counts",colnames(tableOnePyClone))
  colnames(tableOnePyClone)=gsub("major_CN","major_cn",colnames(tableOnePyClone))
  colnames(tableOnePyClone)=gsub("minor_CN","minor_cn",colnames(tableOnePyClone))
  tableAllPyClone=rbind(tableAllPyClone,tableOnePyClone)
}

#run sciClone
library(sciClone)

sc = try(sciClone(vafs=AllelicCountsTableListSciClone,
                  copyNumberCalls=cnvTableListSciClone,
                  sampleNames=sampleNames,
                  minimumDepth=50,
                  #cnCallsAreLog2=TRUE,
                  doClusteringAlongMargins=FALSE,
                  maximumClusters=10
))
if (is.null(sc)) {
  next;
}

outDir="/scratch/cqs/zhaos/workSync/Pierre/20200730_UMI_MuTect2PON_simpleProcessing/sciCloneAllelicCounts/"
outputName=paste(sampleNames,collapse=".")

writeClusterTable(sc, paste0(outDir,outputName,".sciClone.txt"))
sc.plot1d(sc,paste0(outDir,outputName,".sciClone.clusterVis1D.pdf"))
sc.plot2d(sc,paste0(outDir,outputName,".sciClone.clusterVis2D.pdf"))


#make file to run PyClone-VI
tableAllPyClone$major_cn=round(tableAllPyClone$major_cn,0)
tableAllPyClone$minor_cn=round(tableAllPyClone$minor_cn,0)
#set NA to 2,0
tableAllPyClone$major_cn[which(is.na(tableAllPyClone$major_cn))]=2
tableAllPyClone$minor_cn[which(is.na(tableAllPyClone$minor_cn))]=0
write.table(data.frame(tableAllPyClone), paste0(outDir,outputName,".pyCloneVI.input.txt"),sep="\t",quote=FALSE,row.names=FALSE)


############################
#test run pyclone
############################
#source activate pyclone-vi
#export PYTHONPATH=""
#pyclone-vi --help
#cd /scratch/cqs/zhaos/workSync/Pierre/20200730_UMI_MuTect2PON_simpleProcessing/sciCloneAllelicCounts/
#pyclone-vi fit -i Tumor_10658_A.Tumor_10658_B.Tumor_10658_C.pyCloneVI.input.txt -o Tumor_10658_A.Tumor_10658_B.Tumor_10658_C.pyCloneVI.input.txt.pyCloneVI.hd5
#pyclone-vi write-results-file -i Tumor_10658_A.Tumor_10658_B.Tumor_10658_C.pyCloneVI.input.txt.pyCloneVI.hd5 -o Tumor_10658_A.Tumor_10658_B.Tumor_10658_C.pyCloneVI.input.txt.pyCloneVI.txt





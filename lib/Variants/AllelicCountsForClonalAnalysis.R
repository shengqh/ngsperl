

library(data.table)
library(GenomicRanges)

familyInfoFile=parFile1
familyInfoTable=read.delim(familyInfoFile,row.names=1,header=T,as.is=T,check.names=F)

if (!exists("patientFeature")) {
  patientFeature=colnames(familyInfoTable)[2]
}
allPatients=as.character(unique(familyInfoTable[,patientFeature]))

if (!exists("reportingPatientsWithSingleSample")) {
  reportingPatientsWithSingleSample=FALSE
}

#mafFile to define mutation sites for all sample in clonal analysis
#mafFile="/scratch/cqs/zhaos/workSync/Pierre/20200730_UMI_MuTect2PON_simpleProcessing/UMIMutect2PON.filter.allSamples.maf"
mafFile=parFile2
mafTable=fread(mafFile)

#cnvFile="/scratch/cqs/zhaos/Pierre/WES/20200729_UMI_somaticCNV_simpleProcessing/somaticCNV_03_summary/20200729_UMI_somaticCNV_simpleProcessing.allCNV.filter.seg"
cnvFile=parFile3
cnvTable=fread(cnvFile)


AllelicCountsFiles=parSampleFile1
AllelicCountsSampleToFile=read.delim(AllelicCountsFiles,header=FALSE,as.is=T,check.names=F)



sciCloneInputList=data.frame()
pyCloneVIInputList=data.frame()
for (patientOne in allPatients) {
  samplesInSamePatient=row.names(familyInfoTable)[which(familyInfoTable[,patientFeature]==patientOne)]
  
  if (length(samplesInSamePatient)<=1 && !reportingPatientsWithSingleSample) { #skipping patients with only one sample
    print(paste0("Skipping Patient because of only having one sample: ",patientOne))
    next;
  }
  
  print(paste0("Exporting results in Patient: ",patientOne))
  #patient level data
  #samplesInSamePatient=c("Brushing_10658","Tumor_10658_A","Tumor_10658_B","Tumor_10658_C")
  #samplesInSamePatient=c("Tumor_10658_A","Tumor_10658_B","Tumor_10658_C")
  
  ##################################################################################
  #Export positions in samples of this patient to extract AllelicCounts
  ##################################################################################
  mafTablePatientOne=mafTable[Tumor_Sample_Barcode %in% samplesInSamePatient & Variant_Type=="SNP",]
  #posCol=c("Chromosome","Start_Position","End_Position")
  #mafTablePatientOne[,..posCol]
  #mafTablePatientOneGRange=makeGRangesFromDataFrame(mafTablePatientOne[,..posCol],seqnames.field=posCol[1],start.field=posCol[2],end.field=posCol[3],ignore.strand=TRUE)
  #mafTablePatientOneGRange=reduce(mafTablePatientOneGRange)
  posCol=c("Chromosome","Start_Position")
  mafTablePatientOneRangeKey=mafTablePatientOne[,..posCol]
  colnames(mafTablePatientOneRangeKey)=c("CONTIG","POSITION")
  mafTablePatientOneRangeKey=unique(mafTablePatientOneRangeKey)
  mafTablePatientOneRangeKey=mafTablePatientOneRangeKey[order(mafTablePatientOneRangeKey),]
  
  AllelicCountsTableListSciClone=list()
  tableAllPyClone=NULL
  cnvTableListSciClone=list()
  for (sampleOne in samplesInSamePatient) {
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
  
  ########################################################
  #make input file to run PyClone-VI
  ########################################################
  tableAllPyClone$major_cn=round(tableAllPyClone$major_cn,0)
  tableAllPyClone$minor_cn=round(tableAllPyClone$minor_cn,0)
  #set NA to 2,0
  tableAllPyClone$major_cn[which(is.na(tableAllPyClone$major_cn))]=2
  tableAllPyClone$minor_cn[which(is.na(tableAllPyClone$minor_cn))]=0
  pyCloneVIInputFile=paste0(patientOne,".pyCloneVI.input.txt")
  write.table(data.frame(tableAllPyClone), pyCloneVIInputFile,sep="\t",quote=FALSE,row.names=FALSE)
  pyCloneVIInputList=rbind(pyCloneVIInputList,c(patientOne,paste0(getwd(),"/",pyCloneVIInputFile)))
  
  ########################################################
  #make input obj to run sciCLone
  ########################################################
  sciCloneInputFile=paste0(patientOne,".sciClone.input.rds")
  saveRDS(list(AllelicCountsTableListSciClone,cnvTableListSciClone,samplesInSamePatient),sciCloneInputFile)
  sciCloneInputList=rbind(sciCloneInputList,c(patientOne,paste0(getwd(),"/",sciCloneInputFile)))
}


#Save a list of all files as output of this R
colnames(sciCloneInputList)=c("Patient","File")
write.table(sciCloneInputList,paste0(outFile,".sciCloneInputList.txt"),quote=FALSE,row.names = FALSE,sep="\t")

colnames(pyCloneVIInputList)=c("Patient","File")
write.table(pyCloneVIInputList,paste0(outFile,".pyCloneVIInputList.txt"),quote=FALSE,row.names = FALSE,sep="\t")





library(ape)
library(maftools)
library(dplyr)
library(tidyr)

vIdCols=c("Hugo_Symbol", "Protein_Change" ,"Chromosome" ,       "Start_Position"   ,
          "End_Position"   ,   "Reference_Allele" , "Tumor_Seq_Allele2")




familyInfoFile=parFile1
familyInfoTable=read.delim(familyInfoFile,row.names=1,header=T,as.is=T,check.names=F)

if (!exists("patientFeature")) {
  patientFeature=colnames(familyInfoTable)[2]
}
allPatients=as.character(unique(familyInfoTable[,patientFeature]))


#mafFile to define mutation sites for all sample in clonal analysis
#mafFile="/scratch/cqs/zhaos/workSync/Pierre/20200730_UMI_MuTect2PON_simpleProcessing/UMIMutect2PON.filter.allSamples.maf"
mafFile=parFile2
library(data.table)
mafTable=fread(mafFile)



for (patientOne in allPatients) {
  samplesInSamePatient=row.names(familyInfoTable)[which(familyInfoTable[,patientFeature]==patientOne)]
  apeTreeOutputName=paste0(patientOne,".ApePhylogeneticTree")
  
  if (length(samplesInSamePatient)<=2) { #skipping patients with only one or two samples
    print(paste0("Skipping Patient because of only having one or two samples: ",patientOne))
    next;
  }
  
  print(paste0("Exporting results in Patient: ",patientOne))
  #patient level data
  #samplesInSamePatient=c("Brushing_10658","Tumor_10658_A","Tumor_10658_B","Tumor_10658_C")
  #samplesInSamePatient=c("Tumor_10658_A","Tumor_10658_B","Tumor_10658_C")
  
  mafTableOnePatient=mafTable[Tumor_Sample_Barcode %in% samplesInSamePatient]
  variantId=apply(mafTableOnePatient[,..vIdCols],1,function(x) paste(gsub(" ","",x),collapse="_"))
  
  dataForApe=data.frame(Sample=as.character(mafTableOnePatient[["Tumor_Sample_Barcode"]]),Id=variantId)
  dataForApe=dataForApe %>% count(Sample, Id)%>% pivot_wider(names_from=c("Id"),values_from=c("n"),values_fill=0)
  dataForApe=as.data.frame(dataForApe)
  row.names(dataForApe)=dataForApe$Sample
  dataForApe=dataForApe[,-1]
#  dataForApe[1:3,1:3]
  
  
  dataForApeDist=dist.gene(dataForApe)
  #simple neighbor-joining tree 
  dataForApeDistNjTree <- nj(dataForApeDist)
  
  #plot.phylo 
  pdf(paste0(apeTreeOutputName,".pdf"))
  plot(dataForApeDistNjTree ,"unrooted",lab4ut = "axial")
  dev.off()
}

save.image(paste0(outFile,".RData"))


# #Test running ape for ploygentic tree
# 
# library(ape)
# library(maftools)
# library(dplyr)
# library(tidyr)
# 
# #mafFile="D:\\workSync\\Pierre\\20200730_UMI_MuTect2PON_simpleProcessing\\UMIMutect2PON.filter.allSamples.maf"
# #mafObj=read.maf(mafFile,sep = "\t")
# 
# mafReportImage="D:\\workSync\\Pierre\\20200730_UMI_MuTect2PON_simpleProcessing\\UMIMutect2PON.filter.allSamples.maf.report.html.RData"
# load(mafReportImage)
# mafObj=dataForReport[["maf"]]
# 
# 
# samplesInSamePatients=c("Brushing_10658","Tumor_10658_A","Tumor_10658_B","Tumor_10658_C")
# 
# #extract both non-syn and syn variants in these samples for ape
# mafObjSub=subsetMaf(mafObj,tsb =samplesInSamePatients)
# 
# mafTableOnePatient=rbind(mafObjSub@data,mafObjSub@maf.silent)
# 
# if ("vIdCol" %in% names(dataForReport)) {
#   vIdCols=dataForReport$vIdCol
# } else {
#   vIdCols=c("Hugo_Symbol", "Protein_Change" ,"Chromosome" ,       "Start_Position"   ,
#             "End_Position"   ,   "Reference_Allele" , "Tumor_Seq_Allele2")
# }
# variantId=apply(mafTableOnePatient[,..vIdCols],1,function(x) paste(gsub(" ","",x),collapse="_"))
# 
# dataForApe=data.frame(Sample=as.character(mafTableOnePatient[["Tumor_Sample_Barcode"]]),Id=variantId)
# dataForApe=dataForApe %>% count(Sample, Id)%>% pivot_wider(names_from=c("Id"),values_from=c("n"),values_fill=0)
# dataForApe=as.data.frame(dataForApe)
# row.names(dataForApe)=dataForApe$Sample
# dataForApe=dataForApe[,-1]
# dataForApe[1:3,1:3]
# 
# 
# dataForApeDist=dist.gene(dataForApe)
# 
# #simple neighbor-joining tree 
# library(ape)
# 
# dataForApeDistNjTree <- nj(dataForApeDist)
# #plot.phylo 
# plot(dataForApeDistNjTree ,"unrooted")



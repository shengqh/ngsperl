library(data.table)
library(maftools)

args = commandArgs(trailingOnly=TRUE)

inputFile = args[1]
outputFile = args[2]
refBuild = args[3]

annovarFinalTable=fread(sep="\t",cmd=paste0("grep -v '^#' ", inputFile))
tempFile<-paste0(outputFile, ".tmp")
annovarFields = c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene")
tempData<-annovarFinalTable[, colnames(annovarFinalTable) %in% annovarFields, with=FALSE]
leftData<-annovarFinalTable[, !(colnames(annovarFinalTable) %in% annovarFields), with=FALSE]
leftData$uid<-paste(tempData$Chr, tempData$Start, tempData$End, tempData$Ref, tempData$Alt, sep = "_")

tempData$ExonicFunc.refGene = ifelse(test = tempData$Func.refGene == "ncRNA_splicing", yes = "Splice_Site", no = tempData$ExonicFunc.refGene)

write.table(tempData,tempFile,sep="\t",row.names=FALSE, quote = F)
testMafTable=annovarToMaf(tempFile, refBuild=refBuild)
testMafTable$Tumor_Sample_Barcode = ""
class(testMafTable$AAChange)="character"
aaChange<-testMafTable$AAChange
aaChange[aaChange=="character(0)"]<-""
testMafTable$AAChange<-aaChange

testMafTable$uid<-paste(testMafTable$Chromosome, testMafTable$Start_Position, testMafTable$End_Position, testMafTable$Reference_Allele, testMafTable$Tumor_Seq_Allele2, sep = "_")

mafResult = merge(testMafTable, leftData, by = "uid")
mafResult = mafResult[, uid:=NULL]

avsnpIndex = which(grepl("avsnp", colnames(mafResult)))
if (avsnpIndex > 0){
  mafResult$dbSNP_RS = mafResult[,..avsnpIndex]
  mafResult = mafResult[, -avsnpIndex, with=FALSE]
}

write.table(mafResult,outputFile,sep="\t",quote=FALSE,row.names=FALSE)
file.remove(tempFile)


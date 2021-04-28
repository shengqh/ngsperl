

#traform maf file to intervals file for GATK
#https://gatk.broadinstitute.org/hc/en-us/articles/360035531852

mafFileTable=read.delim(parSampleFile1,as.is=T,check.names=F)


#use the largest freq0* file as mafFile to output variant locations
#mafFile="/scratch/cqs/zhaos/Pierre/WES/20200730_UMI_MuTect2PON_simpleProcessing/files_muTect2_02_call_mergeAndMafreport/result/UMIMutect2PON.filter.allSamples.maf"
mafFile=mafFileTable[which.max(file.size(mafFileTable[,1])),1]


library(data.table)

mafTable=fread(mafFile)

#extract SNP only, and unique regions (as samples may have variants in the same regions)
colToExtract=c("Chromosome","Start_Position","End_Position")
mafToIntervals=(mafTable[Variant_Type=="SNP"][,..colToExtract])
mafToIntervals=unique(mafToIntervals)
mafToIntervals=mafToIntervals[order(mafToIntervals),]

mafToIntervals=paste0(mafToIntervals[[1]],":",mafToIntervals[[2]],"-",mafToIntervals[[3]])

#write.table(mafToIntervals,paste0("/scratch/cqs/zhaos/temp/UMIMutect2PON.filter.allSamples.maf.intervals"),row.names=FALSE,quote=FALSE,sep=" ",col.names=FALSE)
#writeLines(mafToIntervals,paste0("/scratch/cqs/zhaos/temp/UMIMutect2PON.filter.allSamples.maf.intervals"))
writeLines(mafToIntervals,paste0(outFile,".maf.intervals"))


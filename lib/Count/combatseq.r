rm(list=ls())
outFile='Mouse_IRI_cardiac_EVs'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2='/nobackup/shah_lab/shengq2/Emeli_Mouse_IRI_cardiac_EVs/20251219_ECa6_ECa9_mm10_SMARTerPicov2//covariance.txt'
parFile3=''


setwd('/nobackup/shah_lab/shengq2/Emeli_Mouse_IRI_cardiac_EVs/20251219_ECa6_ECa9_mm10_SMARTerPicov2/genetable_combatseq/result')

### Parameter setting end ###

library(data.table)
library(dplyr)
library(sva)

count_file_tbl=fread(parSampleFile1, header=FALSE, data.table=FALSE)
countfiles=count_file_tbl$V1

if(parFile2 == ''){
  stop("covariance file not specified")
}

covariate_df=fread(parFile2, header=TRUE, data.table=FALSE)
if(!"batch" %in% colnames(covariate_df)){
  stop(paste0("Covariance file must contain a 'batch' column: ", parFile2))
}

countfile_index = 1
for(countfile_index in c(1:length(countfiles))){
  countfile = countfiles[countfile_index]

  cat("CombatSeq:", countfile, "\n")

  data<-fread(countfile, data.table=FALSE)

  # find the first sample column
  colClass<-sapply(data, class)
  countNotNumIndex<-which((colClass!="numeric" & colClass!="integer" & colnames(data) != "Feature_length" & colnames(data) != "Feature_chr") | grepl("Gene_Id", colnames(data)))
  if (length(countNotNumIndex)==0) {
    index<-1;
    indecies<-c()
  } else {
    index<-max(countNotNumIndex)+1
    indecies<-c(1:(index-1))
  }

  # extract count data
  countData<-data[,c(index:ncol(data)),drop=FALSE]
  countData[is.na(countData)] <- 0
  countData<-round(countData)

  # extract metadata
  metaData<-data[,indecies,drop=FALSE]

  # make sure the covariate_df and countData have matching samples
  stopifnot(nrow(covariate_df) == ncol(countData))
  stopifnot(all(sort(covariate_df$Sample) == sort(colnames(countData))))
  countData=countData[,covariate_df$Sample,drop=FALSE]
  
  # combatseq
  combatData <- ComBat_seq( counts=as.matrix(countData), 
                            batch=as.character(covariate_df$batch), 
                            group=NULL, 
                            covar_mod=NULL)
  
  combat_file=ifelse(grepl("proteincoding.count$", countfile), paste0(outFile, ".combatseq.proteincoding.count"), paste0(outFile, ".combatseq.count"))

  meta_count_data=cbind(metaData, combatData)
  write.table(meta_count_data, file=combat_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

#export session information
sva_version<-paste0("sva,v", packageVersion("sva"))
writeLines(sva_version, paste0(outFile,".sva.version"))

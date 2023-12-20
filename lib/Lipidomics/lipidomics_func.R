library(data.table)
library(dplyr)

read_lipid_tsv_and_filter<-function(filename, filePrefix, qc_name="QC", blank_name="Blank", other_samples=c("AIS", "dBlank")){
  modePrefix=".*[POS|NEG]_"

  header<-data.frame(t(fread(filename, nrows=5, data.table=FALSE)))

  nameCol=min(which("Metabolite name"==header$X5))
  ontologyCol=which("Ontology"==header$X5)
  classCol=which("Class"==header$X1)
  aveCols=which("Average"==header$X4)
  stdevCols=which("Stdev"==header$X4)
  sampleCols=(classCol + 1):(min(aveCols) - 1)
  sampleStartCol=min(sampleCols)

  colnames(header)[1:4]<-header[sampleStartCol-1,c(1:4)]
  
  meta=header[!is.na(header$Class) & !(header$Class %in% c("", "Class")),]
  colnames(meta)[5]="Sample"
  meta$Sample=gsub(modePrefix, "", meta$Sample)
  meta <- meta %>% dplyr::select("Sample", everything())

  write.csv(meta, paste0(filePrefix, ".meta.csv"), row.names=FALSE)

  not_samples=c(qc_name, blank_name, other_samples)
  sample_meta=meta[!(meta[, "Class"] %in% not_samples) & !(meta[, "File type"] %in% not_samples),]
  write.csv(sample_meta, paste0(filePrefix, ".sample_meta.csv"), row.names=FALSE)

  df<-fread(filename, skip=4, data.table=FALSE)
  #convert all integer64 to numeric
  for(i in 1:ncol(df)){
    if(class(df[,i])=="integer64"){
      df[,i]=as.numeric(df[,i])
    }
  }

  colnames(df)[nameCol]<-"Metabolite.name"
  colnames(df)[sampleCols]<-gsub(modePrefix,"",colnames(df)[sampleCols])

  res=df[,c(nameCol, ontologyCol, which(colnames(df) %in% sample_meta$Sample))]

  qc_samples=meta[meta[,"Class"]==qc_name | meta[,"File type"]==qc_name ,"Sample"]
  if(length(qc_samples)==0){
    stop(paste0("No QC sample found for ", qc_name))
  }
  if(length(qc_samples)==1){
    stop(paste("Only 1 QC sample found for", qc_name, ":", qc_samples))
  }

  blank_samples=meta[meta[,"Class"]==blank_name | meta[,"File type"]==blank_name ,"Sample"]
  if(length(blank_samples)==0){
    stop(paste0("No Blank sample found for ", blank_name))
  }
  if(length(blank_samples)==1){
    stop(paste("Only 1 Blank sample found for", blank_name, ":", blank_samples))
  }

  res$QC=apply(df[,qc_samples], 1, mean)
  res$Blank=apply(df[,blank_samples], 1, mean)
  res$Stdev_QC=apply(df[,qc_samples], 1, sd)
  res$RSD_QC=100 * res$Stdev_QC / res$QC
  res$blankQC=res$Blank / res$QC

  write.csv(res, paste0(filePrefix, ".with_QC.csv"), row.names=FALSE)

  filter_tb=data.frame(matrix(ncol=2, nrow=8))
  colnames(filter_tb)=c("Filter", "Count")
  filter_tb[1,]=c("Total features", nrow(res))  

  res=res[res$QC > 0,]
  filter_tb[2,]=c("QC>0", nrow(res))  

  res=res[res$RSD_QC < 25,]
  filter_tb[3,]=c("RSD(QC)<25%", nrow(res))  

  res=res[res$blankQC < 0.1,]
  filter_tb[4,]=c("Blank/QC<0.1", nrow(res))  

  res <- res %>% filter(Metabolite.name != 'Unknown')
  filter_tb[5,]=c('remove "Unknown"', nrow(res))  

  res <- res %>% filter(!grepl('w/o', Metabolite.name))
  filter_tb[6,]=c('remove "w/o MS2"', nrow(res))  

  res <- res %>% filter(Ontology != 'Others') %>% filter(!grepl('RIKEN', Metabolite.name))
  filter_tb[7,]=c('remove other unidentified features', nrow(res))  

  res <- res %>%  filter(!grepl('\\(d7\\)', Metabolite.name)) %>%
                  filter(!grepl('\\-d7\\s', Metabolite.name)) %>%
                  filter(!grepl('\\(d9\\)', Metabolite.name)) %>%
                  filter(!grepl('\\-d9\\s', Metabolite.name)) 

  filter_tb[8,]=c('remove d7/d9', nrow(res)) 

  write.csv(filter_tb, paste0(filePrefix, ".filter.csv"), row.names=FALSE)

  res = res %>% dplyr::select(-c(QC, Blank, Stdev_QC, RSD_QC, blankQC))
  return(list(res=res, all_meta=meta, sample_meta=sample_meta, filter_tb=filter_tb))
}

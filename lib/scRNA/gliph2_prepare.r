rm(list=ls()) 
outFile='AG3669'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/data/h_gelbard_lab/projects/20220508_scRNA_3669/seurat_merge_03_choose_res/result/AG3669.meta.rds'
parFile2='/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_03_convert/result/clonotypes.csv'
parFile3='/data/h_gelbard_lab/projects/20220508_scRNA_3669/arcasHLA_3_merge/result/AG3669.genotypes.tsv'


setwd('/data/h_gelbard_lab/projects/20220508_scRNA_3669/tcr_hla_data/result')

### Parameter setting end ###

library(data.table)
library(reshape2)

meta<-readRDS(parFile1)
cd4cells<-rownames(meta)[grepl("CD4", meta$cell_type)]
cd8cells<-rownames(meta)[grepl("CD8", meta$cell_type)]

subject_tb<-read.table(parSampleFile1, sep='\t', header=F)
samples=unique(subject_tb$V3)

tb<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactor=F)
myoptions<-unlist(split(tb$V1, tb$V2))
gliph2_hla_condition<-unlist(strsplit(myoptions['gliph2_hla_condition'],','))

hla_samples<-subject_tb$V3[subject_tb$V2=='condition' & subject_tb$V1 %in% gliph2_hla_condition]
hla_tb=subject_tb[subject_tb$V3 %in% hla_samples,]

get_subject<-function(sample_names, subject_tb){
  unlist(lapply(sample_names, function(x){ subject_tb$V1[subject_tb$V3 == x & subject_tb$V2 == "Subject"] }))
}

get_condition<-function(sample_names, subject_tb){
  unlist(lapply(sample_names, function(x){ subject_tb$V1[subject_tb$V3 == x & subject_tb$V2 == "condition"] }))
}

get_subject_condition<-function(sample_names, subject_tb){
  subjects = get_subject(sample_names, subject_tb)
  conditions = get_condition(sample_names, subject_tb)
  paste0(subjects, ":", conditions)
}

project_name <- outFile

#load HLA
hla <- fread(parFile3, data.table=F)
hla<-hla[order(hla$subject),]
hla<-hla[hla$subject %in% hla_tb$V3,]
#Modify subject name
hla$subject <- get_subject(hla$subject, hla_tb)

write.table(hla, paste0(project_name, ".hla.txt"), row.names=F, col.names =F, quote=F, sep = "\t")

#load TCR
clonotypes <- fread(parFile2, data.table=F)
x<-clonotypes$cells[1]
cell_types<-lapply(clonotypes$cells, function(x){
  cells<-unlist(strsplit(x, ";"))
  ncd4<-sum(cells %in% cd4cells)
  ncd8<-sum(cells %in% cd8cells)
  ntotal<-sum(cells %in% rownames(meta))
  if(ncd4 >= ntotal * 0.5){
    return(c(ncd4, ncd8, ntotal, "CD4"))
  }
  if(ncd8 >= ntotal * 0.5){
    return(c(ncd4, ncd8, ntotal, "CD8"))
  }
  return(c(ncd4, ncd8, ntotal, "Unknown"))
})

ctdf<-as.data.frame(do.call(rbind, cell_types))
colnames(ctdf)<-c("ncd4", "ncd8", "ntotal", "cell_type")
clonotypes<-cbind(clonotypes, ctdf)
write.csv(clonotypes, paste0(project_name, ".clonotype_celltype.csv"))

#Split CD3
clon_list <- strsplit(clonotypes$cdr3s_aa, ";")
clon_df <- data.frame()
for (i in 1:length(clon_list)){
  clons <- clon_list[[i]]
  if(length(clons) == 2 && grepl("TRA:", clons)){
    clon_df[i,1] <- clons[grep("TRA:", clons)]
    clon_df[i,2] <- clons[grep("TRB:", clons)]
  } else if(length(clons) == 1 && grepl("TRA:", clons)){
    clon_df[i,1] <- clons[grep("TRA:", clons)]
    clon_df[i,2] <- "NA"
  } else if (length(clons) == 1 && grepl("TRB:", clons)){
    clon_df[i,1] <- "NA"
    clon_df[i,2] <- clons[grep("TRB:", clons)]
  } else {
    clon_df[i,1] <- "NA"
    clon_df[i,2] <- "NA"
  }
}

#Data frame
gliph2 <- data.frame(clon_df[,2], clonotypes[, c("v_names","j_names")], clon_df[,1], clonotypes[c(samples, "cell_type")])
colnames(gliph2)[c(1, 4)] <- c("CDR3b", "CDR3a") 

celltype="CD4"
for(celltype in c("CD4", "CD8", "Unknown")){
  ct_df<-gliph2[gliph2$cell_type == celltype, c(1:(ncol(gliph2)-1))]

  #Melt
  gliph2_df <- reshape2::melt(ct_df, id = c("CDR3b", "v_names", "j_names", "CDR3a"), variable.name = "subject", value.name = "count")
  gliph2_df$subject <- get_subject_condition(gliph2_df$subject, subject_tb)

  #Clean
  gliph2_df <- gliph2_df[-c(which(gliph2_df[,"CDR3b"] == "NA")), ]
  gliph2_df[,"CDR3b"] <- gsub("TRB:", "", gliph2_df[,"CDR3b"])
  gliph2_df[,"CDR3a"] <- gsub("TRA:", "", gliph2_df[,"CDR3a"])
  gliph2_df <- gliph2_df[-which(gliph2_df$count == 0), ]

  write.table(gliph2_df, paste0(project_name, ".tcr.", celltype,  ".txt"), row.names=F, col.names =F, quote=F, sep = "\t")
}

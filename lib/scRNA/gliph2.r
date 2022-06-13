outFile="scRNA_5126"
parSampleFile1<-"fileList1.txt"

library(data.table)

project_name <- "s5126"
hla_file <- "/data/h_gelbard_lab/projects/20201228_scRNA_5126/hto_bam_arcasHLA_3_merge/result/scRNA_5126.genotypes.tsv"

#load HLA
hla <- fread(hla_file, data.table=F)

#Modify subject name
hla$subject <- gsub("_.*", "", hla$subject)

write.table(hla, paste(date, project_name, "HLA_gliph2.txt", sep="_"), row.names=F, col.names =F, quote=F, sep = "\t")

library(data.table)
library(reshape)

date <- gsub("-", "", Sys.Date())
project_name <- "s5126"
samples <- c("CTRL_KC", "CTRL_TC", "CTRL_TM", "CTRL_VR", "Pep2_KC", "Pep2_TC", "Pep2_TM", "Pep2_VR")
tcr_file <- "/data/h_gelbard_lab/projects/20201228_scRNA_5126/hto_clonotype_4_convert/result/clonotypes.csv"

#load TCR
clonotypes <- fread(tcr_file, data.table=F)

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
gliph2 <- data.frame(clon_df[,2], clonotypes[, c("TRBV","TRBJ")], clon_df[,1], clonotypes[c(samples)])
colnames(gliph2)[c(1, 4)] <- c("CDR3b", "CDR3a") 

#Melt
gliph2.2 <- melt(gliph2, id = c("CDR3b", "TRBV", "TRBJ", "CDR3a"), 
                 measured = samples)
colnames(gliph2.2)[5:6] <- c("subject", "count") 
gliph2.2[,"subject"] <- gsub("_", ":", gliph2.2[,"subject"])

#Clean
gliph2.2 <- gliph2.2[-c(which(gliph2.2[,"CDR3b"] == "NA")), ]
gliph2.2[,"CDR3b"] <- gsub("TRB:", "", gliph2.2[,"CDR3b"])
gliph2.2[,"CDR3a"] <- gsub("TRA:", "", gliph2.2[,"CDR3a"])
gliph2.2 <- gliph2.2[-which(gliph2.2$count == 0), ]

write.table(gliph2.2, paste(date, project_name, "TCR_gliph2.txt", sep= "_"), row.names=F, col.names =F, quote=F, sep = "\t")


rm(list=ls()) 
outFile='cibersort'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/h_cqs/shengq2/program/collaborations/justin_turner/20210520_CRS/CIBERSORT.R'
parFile2='/nobackup/vickers_lab/projects/20240319_Smith_11160_bulkRNA_human/cibersort/vascular_ref_markers_cpm.rds'
parFile3=''


setwd('/nobackup/vickers_lab/projects/20240319_Smith_11160_bulkRNA_human/cibersort/cibersort/result')

### Parameter setting end ###

##cibersort
library(pheatmap)
library(data.table)
library(edgeR)
library("tools")

source(parFile1)

data_df<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactor=F)

ref_data<-readRDS(parFile2)

has_genemap=FALSE
if(parFile3 != ""){
  has_genemap=TRUE
  geneMap<-read.table(parFile3, sep="\t", header=T)
  geneMap<-geneMap[geneMap$gene_biotype == "protein_coding",]
  geneMap$gene_id_2<-gsub("\\.\\d+","", geneMap$gene_id)
  geneSymbolMap<-split(geneMap$gene_name, geneMap$gene_id_2)
}

idx=1
for(idx in c(1:nrow(data_df))){
  df_name=data_df[idx, 2]
  count_file=data_df[idx, 1]
  cat(df_name, " : ", count_file, "\n")

  if(file_ext(count_file) == "rds"){
    data <- readRDS(count_file)
  }else{
    count <- data.frame(fread(count_file),row.names = 1)
    
    if(has_genemap){
      if ("geneName" %in% colnames(count)){
        library(dplyr)
        geneName=unlist(count$geneName)
        count<-count[,!(colnames(count) %in% c("geneName", "geneDescription"))]
        count2<-aggregate(count, by=list(geneName=geneName), FUN=sum)
        rownames(count2)<-count2$geneName
        count<-count2[,-1]
      }else{
        if(any(rownames(count) %in% geneMap$gene_id_2)){
          #merge gene_symbol
          count<-count[rownames(count) %in% geneMap$gene_id_2,]
          count$geneName<-geneSymbolMap[rownames(count)]
          count<-count[!duplicated(count$geneName),]
          rownames(count)<-count$geneName
          count<-count[,colnames(count) != "geneName"]
        }
      }
    }

    #normalization by CPM
    dge <- DGEList(counts = count,remove.zeros = T)
    tmm_counts <- calcNormFactors(dge)
    #drop <- which(apply(cpm(tmm_counts), 1, max) < 1)
    CPM <- cpm(tmm_counts,prior.count=10)
    data <- CPM
    
    sdata<-cbind(data.frame(Gene=rownames(data)), data)
    write.table(sdata, file=paste0(df_name, ".txt"), sep="\t", row.names=F, quote=F)
  }

  height=max(10, ceiling(ncol(data) / 5))
  
  cat(df_name, " : running CIBERSORT absolute\n")
  cp <- CIBERSORT(ref_data,data,perm = 1000,absolute = T,QN=F)##absolute
  write.csv(cp, file=paste0(df_name, "_CIBERSORT_absolute.cpm.csv"))

  pvalue_index<-which(colnames(cp) == "P-value")
  cpv<-cp[,c(1:(pvalue_index-1))]
  pdf(paste0(df_name, "_CIBERSORT_absolute.cpm.pdf"), width=10, height=height)
  pheatmap(cpv,cluster_cols = F,cluster_rows = F)
  dev.off()
  
  cat(df_name, " : running CIBERSORT relative\n")
  cp <- CIBERSORT(ref_data,data,perm = 1000,absolute = F,QN=F)##relative
  write.csv(cp, file=paste0(df_name,"_CIBERSORT_relative.cpm.csv"))

  pvalue_index<-which(colnames(cp) == "P-value")
  cpv<-cp[,c(1:(pvalue_index-1))]
  pdf(paste0(df_name, "_CIBERSORT_relative.cpm.pdf"), width=10, height=height)
  pheatmap(cpv,cluster_cols = F,cluster_rows = F)
  dev.off()
}

##cibersort
library(pheatmap)
library(data.table)
library(edgeR)
library("tools")

source(parFile1)

data_df<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactor=F)

ref_data<-readRDS(parFile2)

geneMap<-read.table(parFile3, sep="\t", header=T)
geneMap<-geneMap[geneMap$gene_biotype == "protein_coding",]
geneMap$gene_id_2<-gsub("\\.\\d+","", geneMap$gene_id)
geneSymbolMap<-split(geneMap$gene_name, geneMap$gene_id_2)

idx=1
for(idx in c(1:nrow(data_df))){
  df_name=data_df[idx, 2]
  count_file=data_df[idx, 1]
  cat(df_name, " : ", count_file, "\n")

  if(file_ext(count_file) == "csv"){
    count <- read.csv(count_file,row.names = 1)
  }else{
    count <- read.table(count_file,row.names = 1, sep="\t", header=T)
  }
  
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
  
  #normalization by CPM
  dge <- DGEList(counts = count,remove.zeros = T)
  tmm_counts <- calcNormFactors(dge)
  #drop <- which(apply(cpm(tmm_counts), 1, max) < 1)
  CPM <- cpm(tmm_counts,prior.count=10)
  data <- CPM
  
  sdata<-cbind(data.frame(Gene=rownames(data)), data)
  write.table(sdata, file=paste0(df_name, ".txt"), sep="\t", row.names=F, quote=F)
  
  height=max(10, ceiling(ncol(data) / 5))
  
  cp <- CIBERSORT(ref_data,data,perm = 1000,absolute = T,QN=F)##absolute
  write.csv(cp, file=paste0(df_name, "_CIBERSORTx_absolute.cpm.csv"))

  pvalue_index<-which(colnames(cp) == "P-value")
  cpv<-cp[,c(1:(pvalue_index-1))]
  pdf(paste0(df_name, "_CIBERSORTx_absolute.cpm.pdf"), width=10, height=height)
  pheatmap(cpv,cluster_cols = F,cluster_rows = F)
  dev.off()
  
  cp <- CIBERSORT(ref_data,data,perm = 1000,absolute = F,QN=F)##relative
  write.csv(cp, file=paste0(df_name,"_CIBERSORTx_relative.cpm.csv"))

  pvalue_index<-which(colnames(cp) == "P-value")
  cpv<-cp[,c(1:(pvalue_index-1))]
  pdf(paste0(df_name, "_CIBERSORTx_relative.cpm.pdf"), width=10, height=height)
  pheatmap(cpv,cluster_cols = F,cluster_rows = F)
  dev.off()
}

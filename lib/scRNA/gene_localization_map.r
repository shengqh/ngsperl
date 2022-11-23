rm(list=ls()) 
outFile='mouse_8870'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/scratch/jbrown_lab/shengq2/projects/20221117_scRNA_8870_mouse/seurat_sct_merge/result/mouse_8870.final.rds'
parFile2=''
parFile3=''


setwd('/scratch/jbrown_lab/shengq2/projects/20221117_scRNA_8870_mouse/seurat_sct_merge_gene_localization_map/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(ggpubr)

obj<-read_object(parFile1)

if(file.exists(parSampleFile2)){
  groups_tbl<-read.table(parSampleFile2, sep="\t", stringsAsFactors = F)
  groups=split(groups_tbl$V2, groups_tbl$V1)
  obj$group = unlist(groups[obj$orig.ident])

  ngroup=length(unique(groups_tbl$V2))
  has_group = TRUE
}else{
  obj$group = "all"
  has_group = FALSE
}

#using RNA assay for visualization
DefaultAssay(obj)<-"RNA"

filelist=NULL
genes_tbl<-read.table(parSampleFile1, sep="\t", stringsAsFactors = F)
genes<-genes_tbl$V1
genes<-gsub('\\s+','',genes)
miss_genes<-genes[!(genes %in% rownames(obj))]

if(length(miss_genes) > 0){
  writeLines(miss_genes, con="miss_genes.txt")
  warning(paste0("There is missing genes:", paste0(miss_genes, ",")))
  genes<-genes[genes %in% rownames(obj),]
}
  
gene=genes_tbl$V1[1]
for (gene in genes_tbl$V1){
  g1<-FeaturePlot(obj, gene, cols=c("lightgray", "red"), order=TRUE, pt.size=1)
  
  pngfile = paste0(outFile, ".", gene, ".all.png")
  png(filename=pngfile, width=1300, height=1100, res=300)
  print(g1)
  dev.off()
  filelist=rbind(filelist, data.frame("file"=paste0(getwd(), "/", pngfile), "gene"=gene, "type"="all"))

  if(has_group){
    g2<-FeaturePlot(obj, gene, cols=c("lightgray", "red"), order=TRUE, pt.size=1, split.by="group")
    pngfile = paste0(outFile, ".", gene, ".group.png")
    png(filename=pngfile, width= ngroup * 1200, height=1200, res=300)
    print(g2)
    dev.off()
    filelist=rbind(filelist, data.frame("file"=paste0(getwd(), "/", pngfile), "gene"=gene, "type"="group"))
  }
}

write.csv(filelist, paste0(outFile, ".figure.files.csv"), row.names=F)

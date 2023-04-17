rm(list=ls()) 
outFile='mouse_8363'
parSampleFile1=''
parSampleFile2=''
parSampleFile3=''
parFile1='/scratch/jbrown_lab/shengq2/projects/20220630_scRNA_8363_mouse/seurat_sct_merge_multires_03_choose/result/mouse_8363.final.rds'
parFile2=''
parFile3=''
genes='Cd68;Cd163;Adgre1;Ccr2;Ccr5;Cx3cr1;Brca1;Rad51;Trl9;Tmem173;cgas;Mx1;Ifit1;Irf7;'; celltype_name='cell_type'; cluster_name='seurat_cell_type'; dotPlotOnly=0;

setwd('/scratch/jbrown_lab/shengq2/projects/20220630_scRNA_8363_mouse/seurat_sct_merge_multires_03_choose_genes/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(ggpubr)

genes<-unlist(strsplit(genes, ";"))
genes<-unique(genes)

obj<-read_object(parFile1)

DefaultAssay(obj)<-"RNA"

getCells<-function(curassay, assay_name, genes){
  inrawcount<-genes %in% rownames(curassay)
  cells<-lapply(c(1:length(genes)), function(x){
    inraw<-inrawcount[x]
    gene<-genes[x]
    #cat(inraw, gene, "\n")
    if (!inraw){
      return(c(0, ncol(curassay)))
    }else{
      incell<-sum(curassay[gene,] > 0)
      notincell<-ncol(curassay)-incell
      return(c(incell, notincell))
    }
  })
  cells<-data.frame(matrix(unlist(cells), nrow=length(cells), byrow=T))
  result<-data.frame(X1 =inrawcount, X2=cells$X1, X3=cells$X2)
  colnames(result)<-c(paste0("InAssay_",assay_name), paste0("InCell_", assay_name), paste0("NotInCell_", assay_name))
  rownames(result)<-genes
  return(result)
}

rawgenes<-getCells(obj@assays$RNA, "RNA", genes)
activegenes<-getCells(obj@assays[[obj@active.assay]], obj@active.assay, genes)

geneinfo<-cbind(rawgenes, activegenes)
write.csv(geneinfo, file="gene_summary.csv")

genes<-genes[genes %in% rownames(obj)]

if(parFile3 != ""){
  clusterDf<-read.csv(parFile3, stringsAsFactors = F)
}else{
  clusterDf<-obj@meta.data
}
clusters<-clusterDf[,cluster_name]

caCount<-table(clusters)
clusterDf$caCount<-caCount[clusters]

clusterDf<-clusterDf[order(clusterDf$caCount, clusterDf[,cluster_name], decreasing = c(1,0)),]

clusters<-factor(clusters, levels=unique(clusterDf[,cluster_name]))

obj$final_seurat_clusters<-clusters

samples<-unique(obj$orig.ident)

if(!exists("dotPlotOnly")){
  dotPlotOnly<-FALSE
}

if(!dotPlotOnly){
  gene=genes[1]
  ncol=ceiling(sqrt(1 + length(samples)))
  nrow=ceiling((1 + length(samples)) / ncol)
  for(gene in genes){
    p1<-MyDimPlot(obj, reduction = "umap", label=T, group.by="final_seurat_clusters") + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5))
    lst<-list("Cluster" = p1)
    sample<-samples[1]
    for (sample in samples){
      sobj<-subset(obj, cells=colnames(obj)[obj$orig.ident==sample])
      pgene<-MyFeaturePlot(object = sobj, features=gene)  + NoLegend() + ggtitle(sample) + theme(plot.title = element_text(hjust=0.5))
      lst[[sample]]=pgene
    }
    
    png(filename=paste0(outFile, ".", gene, ".png"), width= ncol * 3000, height=nrow * 3000, res=300)
    g<-ggarrange(plotlist=lst, ncol=ncol, nrow=nrow)
    print(g)
    dev.off()
  }
}

pdf(file=paste0(outFile, ".dot.pdf"), width=max(length(genes) * 0.4, 10), height=10, onefile = T)

alltitle = ifelse(length(samples) == 1, "", "All samples")

#png(filename=paste0(outFile, ".dot.png"), width=max(length(genes) * 100, 5000), height=2500, res=300)
p<-DotPlot(obj, assay = "RNA", group.by="final_seurat_clusters", features=genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
  xlab("genes") + ggtitle(alltitle) + theme(plot.title = element_text(hjust = 0.5))
print(p)

if (length(samples) > 1) {
  for (sample in samples){
    sobj<-subset(obj, cells=colnames(obj)[obj$orig.ident==sample])
    p<-DotPlot(sobj, assay = "RNA", group.by="final_seurat_clusters", features=genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
      xlab("genes") + ggtitle(paste0("Sample ", sample)) + theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
}
dev.off()

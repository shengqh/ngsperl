rm(list=ls()) 
outFile='SADIE_pbmc'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/shah_lab/shengq2/20240304_mona_scRNA_SADIE/data/pbmc.v4.DE.rds'
parFile2='/nobackup/shah_lab/shengq2/20240304_mona_scRNA_SADIE/20240305_DE_fold1.2_pbmc/files_edgeR_inCluster_bySample/result/SADIE_pbmc.edgeR.files.csv'
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20240304_mona_scRNA_SADIE/20240305_DE_fold1.2_pbmc/files_edgeR_inCluster_bySample_vis/result')

### Parameter setting end ###

source("scRNA_func.r")
library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(R.utils)
library(reshape2)
library(patchwork)

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
bBetweenCluster<-is_one(myoptions$bBetweenCluster)
cluster_name=myoptions$cluster_name
DE_by_cell=is_one(myoptions$DE_by_cell)
reduction=myoptions$reduction
sample_column=myoptions$sample_column

if(!exists("obj")){
  obj<-read_object(parFile1, parFile3, cluster_name)
  if(cluster_name == "bulk"){
    obj@meta.data$bulk = "bulk"
  }else{
    obj@meta.data[,cluster_name]=gsub("^\\s+", "",obj@meta.data[,cluster_name])
  }
  if(sample_column != "orig.ident"){
    obj@meta.data$orig.ident = obj@meta.data[,sample_column]
  }
}
clusterDf<-obj@meta.data

edgeRres<-read.csv(parFile2, stringsAsFactors = F, row.names=1)
edgeRres$prefix<-basename(edgeRres$prefix)
edgeRfolder<-dirname(parFile2)
rownames(edgeRres)<-edgeRres$prefix

if(!bBetweenCluster){
  if(all(is.null(levels(obj@meta.data[, cluster_name])))){
    if("seurat_clusters" %in% colnames(obj@meta.data)){
      if(length(unique(obj$seurat_clusters)) == length(unique(obj@meta.data[, cluster_name]))){
        df<-data.frame(c1=obj$seurat_clusters, c2=as.character(unlist(obj@meta.data[, cluster_name])))
        df<-unique(df)
        df<-df[order(df$c1),]
        obj@meta.data[, cluster_name]<-factor(as.character(unlist(obj@meta.data[, cluster_name])), levels=unique(df$c2))
      }else{
        obj@meta.data[, cluster_name]<-factor_by_count(obj@meta.data[, cluster_name])
      }
    }else{
      obj@meta.data[, cluster_name]<-factor_by_count(obj@meta.data[, cluster_name])
    }
  }
}

all_sigout<-NULL

result<-NULL
prefix<-rownames(edgeRres)[2]
for (prefix in rownames(edgeRres)){
  cat("Processing ", prefix, "\n")
  comparison<-edgeRres[prefix, "comparison"]
  sigGenenameFile<-paste0(edgeRfolder, "/", edgeRres[prefix, "sigGenenameFile"])
  cellType<-edgeRres[prefix, "cellType"]
  deFile=paste0(edgeRfolder, "/", edgeRres[prefix, "deFile"])
  totalGene=length(readLines(deFile))-1
  sigGene=length(readLines(sigGenenameFile))
  visFile=""
  if (sigGene > 0){
    visFile=paste0(prefix, ".sig_genename.pdf")
    #if(!file.exists(visFile)){
      siggenes<-read.table(sigGenenameFile, sep="\t", stringsAsFactors = F)
      colnames(siggenes)<-"gene"
      
      sigoutFile<-paste0(edgeRfolder, "/", edgeRres[prefix, "sigFile"])
      sigout<-read.csv(sigoutFile, header=T, stringsAsFactors = F, row.names=1)
      sigout$comparison=prefix
      
      all_sigout<-rbind(all_sigout, sigout)

      designFile<-paste0(edgeRfolder, "/", edgeRres[prefix, "designFile"])
      design_data<-read.csv(designFile, stringsAsFactors = F, header=T)
      
      designUniq<-unique(design_data[,c("Group", "DisplayGroup")])
      rownames(designUniq)<-designUniq$Group
      
      controlGroup<-designUniq["control","DisplayGroup"]
      sampleGroup<-designUniq["sample","DisplayGroup"]

      groupColors<-c("blue", "red")
      names(groupColors)<-c(controlGroup, sampleGroup)

      if(DE_by_cell){
        cell_obj<-subset(obj, cells=design_data$Cell)
        cell_obj$Group=design_data$Group
        cell_obj$DisplayGroup=design_data$DisplayGroup
      }else{
        #pseudo_bulk
        cells<-clusterDf[clusterDf[,cluster_name] == cellType,]
        cells<-cells[cells$orig.ident %in% design_data$Sample,]
        cell_obj<-subset(obj, cells=rownames(cells))

        gmap<-unlist(split(design_data$Group, design_data$Sample))
        gdismap<-unlist(split(design_data$DisplayGroup, design_data$Sample))

        cell_obj@meta.data$Group=gmap[cell_obj$orig.ident]
        cell_obj@meta.data$DisplayGroup=factor(gdismap[cell_obj$orig.ident], levels=names(groupColors))
      
        cpmFile=paste0(edgeRfolder, "/", edgeRres[prefix, "cpmFile"])
        log_cpm<-read.csv(cpmFile, row.names=1)
      }
      
      if("percent.mt" %in% colnames(cell_obj)){
        vars.to.regress="percent.mt"
      }else{
        vars.to.regress=NULL
      }
      cell_obj = ScaleData(cell_obj, features=siggenes$gene, assay="RNA", vars.to.regress=vars.to.regress)
      g<-MyDoHeatMap(cell_obj, assay="RNA", features=siggenes$gene, group.by="DisplayGroup")
      png(paste0(prefix, ".sig_gene.heatmap.png"), 
          width=4000, 
          height=get_heatmap_height(nrow(siggenes)), 
          res=300)
      print(g)
      dev.off()

      coords<-data.frame(Embeddings(cell_obj, reduction=reduction))

      colnames(coords)<-c("UMAP_1", "UMAP_2")

      xlim<-c(min(coords$UMAP_1-0.1), max(coords$UMAP_1+0.1))
      ylim<-c(min(coords$UMAP_2-0.1), max(coords$UMAP_2+0.1))
      
      visFile=paste0(prefix, ".sig_genename.pdf")
      
      if(bBetweenCluster){
        width<-10
      }else{
        nsamples<-length(unique(cell_obj$orig.ident))
        width<-10
        if(nsamples>10){
          width=width+5
        }else if(!DE_by_cell){
          width=width+5
        }
      }

      topN = min(100, nrow(siggenes))
      topNgenes<-siggenes$gene[1:topN]

      pdf(file=visFile, onefile = T, width=width, height=10)
      siggene<-topNgenes[1]
      for (sig_gene in topNgenes){
        p<-get_sig_gene_figure(cell_obj, sigout, design_data, sig_gene, DE_by_cell=DE_by_cell, is_between_cluster=bBetweenCluster, log_cpm=log_cpm)

        print(p)
        #break
      }
      dev.off()
    #}
  }
  curDF<-data.frame("prefix"=prefix, "sigGeneVisFile"=visFile, "sigGene"=sigGene, "totalGene"=totalGene, "cluster"=cellType, "comparison"=comparison)
  if(is.null(result)){
    result<-curDF
  }else{
    result<-rbind(result, curDF)
  }
}

write.csv(result, file=paste0(outFile, ".vis.files.csv"), quote=F)

result$sigRate<-result$sigGene * 100.0 / result$totalGene

if(!bBetweenCluster){
  allcoords<-data.frame(Embeddings(obj, reduction="umap"))
  allcoords$Cluster=obj[[cluster_name]]
  
  comp=unique(result$comparison)[1]
  for (comp in unique(result$comparison)) {
    compRes = result[result$comparison == comp,]
    rateMap=unlist(split(compRes$sigRate, compRes$cluster))
    
    obj@meta.data$sigRate=unlist(rateMap[as.character(unlist(obj[[cluster_name]]))])
    obj@meta.data$sigRate[is.na(obj$sigRate)]<-0
    
    pdf(paste0(outFile, ".", comp, ".sigGenePerc.pdf"), width=14, height=7)
    p1<-get_dim_plot_labelby(obj, label.by = cluster_name) + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5))
    p2<-MyFeaturePlot(obj, feature="sigRate", cols=c("lightgrey", "red")) + ggtitle("Percentage of DE genes in each cluster") + theme(plot.title = element_text(hjust=0.5))
    g<-ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))
    print(g)
    dev.off()
  }
}

write.csv(all_sigout, file=paste0(outFile, ".allsigout.csv"), quote=F)

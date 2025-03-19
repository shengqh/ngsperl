rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose/result/combined.final.rds'
parFile2='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose_edgeR_inCluster_byCell/result/combined.edgeR.files.csv'
parFile3='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose/result/combined.meta.rds'


setwd('/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose_edgeR_inCluster_byCell_vis/result')

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

if(is.null(myoptions$edgeR_suffix)){
  detail_folder = paste0(outFile,  ".edgeR_vis/")
}else{
  detail_folder = paste0(outFile, myoptions$edgeR_suffix, ".vis/")
}
if(!file.exists(detail_folder)){
  dir.create(detail_folder)
}
all_sigout<-NULL

result<-NULL
comp<-rownames(edgeRres)[2]
for (comp in rownames(edgeRres)){
  cat("Processing ", comp, "\n")
  comparison<-edgeRres[comp, "comparison"]
  sigGenenameFile<-paste0(edgeRfolder, "/", edgeRres[comp, "sigGenenameFile"])
  cellType<-edgeRres[comp, "cellType"]
  deFile=paste0(edgeRfolder, "/", edgeRres[comp, "deFile"])
  totalGene=length(readLines(deFile))-1
  sigGene=length(readLines(sigGenenameFile))
  visFile=""

  top_png=""
  prefix=paste0(detail_folder, comp)
  if (sigGene > 0){
    visFile=paste0(prefix, ".sig_genename.pdf")
    #if(!file.exists(visFile)){
      siggenes<-read.table(sigGenenameFile, sep="\t", stringsAsFactors = F)
      colnames(siggenes)<-"gene"
      
      sigoutFile<-paste0(edgeRfolder, "/", edgeRres[comp, "sigFile"])
      sigout<-read.csv(sigoutFile, header=T, stringsAsFactors = F, row.names=1)
      sigout$comparison=prefix
      
      all_sigout<-rbind(all_sigout, sigout)

      designFile<-paste0(edgeRfolder, "/", edgeRres[comp, "designFile"])
      designdata<-read.csv(designFile, stringsAsFactors = F, header=T) |>
        tibble::column_to_rownames("Cell")
      
      groupColors<-get_group_colors_from_designdata(designdata)

      if(DE_by_cell){
        cell_obj<-subset(obj, cells=rownames(designdata))
        # Using cells in subset cannot guarantee the order of cells in cell_obj is the same as in designdata
        # So we must match the order when assign Group and DisplayGroup
        cell_obj$Group=designdata[colnames(cell_obj), "Group"]
        cell_obj$DisplayGroup=designdata[colnames(cell_obj), "DisplayGroup"]
      }else{
        #pseudo_bulk
        cells<-clusterDf[clusterDf[,cluster_name] == cellType,]
        cells<-cells[cells$orig.ident %in% designdata$Sample,]
        cell_obj<-subset(obj, cells=rownames(cells))

        gmap<-unlist(split(designdata$Group, designdata$Sample))
        gdismap<-unlist(split(designdata$DisplayGroup, designdata$Sample))

        cell_obj@meta.data$Group=gmap[cell_obj$orig.ident]
        cell_obj@meta.data$DisplayGroup=factor(gdismap[cell_obj$orig.ident], levels=names(groupColors))
      
        cpmFile=paste0(edgeRfolder, "/", edgeRres[comp, "cpmFile"])
        log_cpm<-read.csv(cpmFile, row.names=1)
      }
      
      if("percent.mt" %in% colnames(cell_obj)){
        vars.to.regress="percent.mt"
      }else{
        vars.to.regress=NULL
      }

      #we have to use more than one genes for scaling
      scale_genes=unique(c(rownames(cell_obj)[1:5], siggenes$gene))
      cell_obj = myScaleData(cell_obj, features=scale_genes, assay="RNA", vars.to.regress=vars.to.regress)

      heatmap_height=min(8, max(4, 0.4*sigGene))
      if(grepl("^5", packageVersion("Seurat"))){
        g<-MyDoHeatMap( cell_obj, 
                        assay="RNA", 
                        features=siggenes$gene, 
                        group.by="DisplayGroup", 
                        group.colors=groupColors, 
                        angle=0,
                        group.bar.height=0.05,
                        vjust=0.1) + guides(color = "none")
      }else{
        g<-MyDoHeatMap( cell_obj, 
                        assay="RNA", 
                        features=siggenes$gene, 
                        group.by="DisplayGroup", 
                        group.colors=groupColors, 
                        angle=0,
                        group.bar.height=0.05) + guides(color = "none")
      }
      if(sigGene > 20){
        g<-g + theme(axis.text.y=element_blank())
      }
      ggsave(paste0(prefix, ".sig_gene.heatmap.png"), g, width=8, height=heatmap_height, dpi=300, units="in", bg="white")

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
      sig_gene<-topNgenes[1]
      is_first=TRUE
      for (sig_gene in topNgenes){
        p<-get_sig_gene_figure( cell_obj, 
                                sigout, 
                                designdata, 
                                sig_gene, 
                                DE_by_cell=DE_by_cell, 
                                is_between_cluster=bBetweenCluster, 
                                log_cpm=log_cpm)

        print(p)
        if(is_first){
          top_png = paste0(prefix, ".", sig_gene, ".png")
          ggsave(top_png, p, width=8, height=6, dpi=300, units="in", bg="white")
          is_first=FALSE
        }
        #break
      }
      dev.off()
    #}
  }
  curDF<-data.frame("prefix"=prefix, "sigGeneVisFile"=visFile, "sigGene"=sigGene, "totalGene"=totalGene, "cluster"=cellType, "comparison"=comparison, "top_png"=top_png)
  if(is.null(result)){
    result<-curDF
  }else{
    result<-rbind(result, curDF)
  }
}

result_csv=paste0(detail_folder, outFile, ".vis.files.csv")
write.csv(result, file=result_csv, quote=F)

if(0){
  result=data.frame(fread(result_csv), row.names=1)
}
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
    
    if(cluster_name == "seurat_cell_type"){
      p1<-get_dim_plot(obj, group.by="seurat_clusters", label.by = cluster_name) 
    }else{
      p1<-get_dim_plot_labelby(obj, label.by = cluster_name)
    }
    p1 <- p1 + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5)) + guides(fill=guide_legend(ncol =1))
    p2<-MyFeaturePlot(obj, feature="sigRate", cols=c("lightgrey", "red"), raster=FALSE) + ggtitle("Percentage of DE genes in each cluster") + theme(plot.title = element_text(hjust=0.5))
    g<-ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))
    ggsave(paste0(detail_folder, outFile, ".", comp, ".sigGenePerc.png"), g, width=14, height=7, dpi=300, units="in", bg="white")
  }
}

write.csv(all_sigout, file=paste0(detail_folder, outFile, ".allsigout.csv"), quote=F)

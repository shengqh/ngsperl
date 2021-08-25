
library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(R.utils)

finalList<-readRDS(parFile1)
obj<-finalList$obj

edgeRres<-read.csv(parFile2, stringsAsFactors = F, row.names=1)
edgeRfolder<-dirname(parFile2)
rownames(edgeRres)<-edgeRres$prefix

clusterDf<-read.csv(parFile3, stringsAsFactors = F, row.names=1)
obj[[cluster_name]]<-clusterDf[names(obj$orig.ident), cluster_name]

df<-data.frame(c1=obj$seurat_clusters, c2=obj[[cluster_name]])
df<-unique(df)
df<-df[order(df$c1),]
obj[[cluster_name]]<-factor(unlist(obj[[cluster_name]]), levels=unique(df[,cluster_name]))

result<-NULL
prefix<-rownames(edgeRres)[2]
for (prefix in rownames(edgeRres)){
  cat("Processing ", prefix, "\n")
  comparison<-edgeRres[prefix, "comparison"]
  sigGenenameFile<-paste0(edgeRfolder, "/", edgeRres[prefix, "sigGenenameFile"])
  cellType<-edgeRres[prefix, "cellType"]
  deFile=gsub(".sig_genename.txt", ".csv", sigGenenameFile)
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
      
      designFile<-paste0(edgeRfolder, "/", edgeRres[prefix, "designFile"])
      design<-read.csv(designFile, stringsAsFactors = F, header=T)
      
      all_cells<-design$Cell
      
      cell_obj<-subset(obj, cells=design$Cell)
      cell_obj$Group=design$Group
      cell_obj$DisplayGroup=design$DisplayGroup
      
      designUniq<-unique(design[,c("Group", "DisplayGroup")])
      rownames(designUniq)<-designUniq$Group
      
      controlGroup<-designUniq["control","DisplayGroup"]
      sampleGroup<-designUniq["sample","DisplayGroup"]
      
      coords<-data.frame(cell_obj@reductions$umap@cell.embeddings)
      xlim<-c(min(coords$UMAP_1-0.1), max(coords$UMAP_1+0.1))
      ylim<-c(min(coords$UMAP_2-0.1), max(coords$UMAP_2+0.1))
      
      visFile=paste0(prefix, ".sig_genename.pdf")
      pdf(file=visFile, onefile = T, width=10, height=10)
      siggene<-siggenes$gene[1]
      for (siggene in siggenes$gene){
        logFC<-sigout[siggene, "logFC"]
        FDR<-sigout[siggene,"FDR"]
        
        geneexp=FetchData(cell_obj,vars=c(siggene))
        colnames(geneexp)<-"Gene"
        colorRange<-c(min(geneexp), max(geneexp))
        fix.sc <- scale_color_gradientn(colors=c("lightgrey", "blue"), limits = colorRange)
        
        geneexp$Group<-cell_obj$DisplayGroup
        geneexp$Sample<-cell_obj$orig.ident
        
        title<-paste0(siggene, ' : logFC = ', round(logFC, 2), ", FDR = ", formatC(FDR, format = "e", digits = 2))
        
        if(bBetweenCluster){
          p0<-ggplot(geneexp, aes(x=Group, y=Gene, col=Group)) + geom_violin() + geom_jitter(width = 0.2) + facet_grid(~Sample) + theme_bw() + NoLegend() + xlab("") + ylab("Gene Expression") + theme(strip.background=element_blank())
          
          p1<-DimPlot(cell_obj, reduction = "umap", label=T, group.by="DisplayGroup") + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5)) + xlim(xlim) + ylim(ylim)
          
          p2<-FeaturePlot(object = cell_obj, features=as.character(siggene), order=T)
        }else{
          p0<-ggplot(geneexp, aes(x="1", y=Gene, color=Group)) + geom_violin() + geom_jitter(width = 0.2) + facet_grid(~Sample) + theme_bw() + xlab("") + ylab("Gene Expression") + theme(strip.background=element_blank(), axis.text.x = element_blank())
          
          subcells<-colnames(cell_obj)[cell_obj$DisplayGroup == controlGroup]
          subobj<-subset(cell_obj, cells=subcells)
          p2<-FeaturePlot(object = subobj, features=siggene, order=T) + ggtitle(paste0("Control: ", controlGroup))
          p2<-suppressMessages(expr = p2 + xlim(xlim) + ylim(ylim) + fix.sc)
          
          subcells<-colnames(cell_obj)[cell_obj$DisplayGroup == sampleGroup]
          subobj<-subset(cell_obj, cells=subcells)
          p1<-FeaturePlot(object = subobj, features=siggene, order=T) + ggtitle(paste0("Sample: ", sampleGroup))
          p1<-suppressMessages(expr = p1  + xlim(xlim) + ylim(ylim) + fix.sc)
        }
        p<-ggarrange(p0,                                                 # First row with scatter plot
                     ggarrange(p1, p2, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
                     nrow = 2, 
                     labels = "A"                                        # Labels of the scatter plot
        ) 
        g<-ggpubr::annotate_figure(
          p = p,
          top = ggpubr::text_grob(label = title, face = 'bold', size=20)
        )
        print(g)
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
  allcoords<-data.frame(obj@reductions$umap@cell.embeddings)
  allcoords$Cluster=obj[[cluster_name]]
  
  for (comp in unique(result$comparison)) {
    compRes = result[result$comparison == comp,]
    rownames(compRes)=compRes$cluster
    
    obj$sigRate=compRes[unlist(obj[[cluster_name]]), "sigRate"]
    
    pdf(paste0(outFile, ".", comp, ".sigGenePerc.pdf"), width=14, height=7)
    p1<-DimPlot(obj, group.by = cluster_name, label=T) + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5))
    p2<-FeaturePlot(obj, feature="sigRate", cols=c("lightgrey", "red")) + ggtitle("Percentage of DE genes in each cluster") + theme(plot.title = element_text(hjust=0.5))
    g<-ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))
    print(g)
    dev.off()
  }
}


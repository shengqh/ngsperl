#rm(list=ls()) 
outFile='AK6383'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose/result/AK6383.final.rds'
parFile2='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose_edgeR_inCluster_bySample/result/AK6383.edgeR.files.csv'
parFile3='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose/result/AK6383.meta.rds'


setwd('C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose_edgeR_inCluster_bySample_vis/result')

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
bBetweenCluster<-ifelse(myoptions$bBetweenCluster == "0", FALSE, TRUE)
cluster_name=myoptions$cluster_name
DE_by_cell=ifelse(myoptions$DE_by_cell == '1', TRUE, FALSE)

if(!exists("obj")){
  obj<-read_object(parFile1, parFile3, cluster_name)
}
clusterDf<-obj@meta.data

edgeRres<-read.csv(parFile2, stringsAsFactors = F, row.names=1)
edgeRfolder<-dirname(parFile2)
rownames(edgeRres)<-edgeRres$prefix

df<-data.frame(c1=obj$seurat_clusters, c2=obj[[cluster_name]])
df<-unique(df)
df<-df[order(df$c1),]
obj[[cluster_name]]<-factor(unlist(obj[[cluster_name]]), levels=unique(df[,cluster_name]))

all_sigout<-NULL

result<-NULL
prefix<-rownames(edgeRres)[4]
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
      sigout$comparison=prefix
      
      all_sigout<-rbind(all_sigout, sigout)

      designFile<-paste0(edgeRfolder, "/", edgeRres[prefix, "designFile"])
      design<-read.csv(designFile, stringsAsFactors = F, header=T)
      
      designUniq<-unique(design[,c("Group", "DisplayGroup")])
      rownames(designUniq)<-designUniq$Group
      
      controlGroup<-designUniq["control","DisplayGroup"]
      sampleGroup<-designUniq["sample","DisplayGroup"]

      groupColors<-c("blue", "red")
      names(groupColors)<-c(controlGroup, sampleGroup)

      if(DE_by_cell){
        cell_obj<-subset(obj, cells=design$Cell)
        cell_obj$Group=design$Group
        cell_obj$DisplayGroup=design$DisplayGroup
      }else{
        #pseudo_bulk
        gmap<-unlist(split(design$Group, design$Sample))
        gdismap<-unlist(split(design$DisplayGroup, design$Sample))
        
        cells<-clusterDf[clusterDf[,cluster_name] == cellType,]
        cells<-cells[cells$orig.ident %in% names(gmap),]
        cell_obj<-subset(obj, cells=rownames(cells))
        
        cell_obj$Group=gmap[cell_obj$orig.ident]
        cell_obj$DisplayGroup=factor(gdismap[cell_obj$orig.ident], levels=names(groupColors))
      
        cpmFile=paste0(edgeRfolder, "/", edgeRres[prefix, "cpmFile"])
        cpmvalues<-read.csv(cpmFile)
        melt_cpm<-melt(cpmvalues, id.vars = 'X')
        colnames(melt_cpm)<-c("Gene", "Sample", "CPM")
        melt_cpm$Group<-factor(gdismap[melt_cpm$Sample], levels=names(groupColors))
      }
      
      coords<-data.frame(cell_obj@reductions$umap@cell.embeddings)
      xlim<-c(min(coords$UMAP_1-0.1), max(coords$UMAP_1+0.1))
      ylim<-c(min(coords$UMAP_2-0.1), max(coords$UMAP_2+0.1))
      
      visFile=paste0(prefix, ".sig_genename.pdf")
      
      nsamples<-length(unique(cell_obj$orig.ident))
      width<-10
      if(nsamples>10){
        width=width+5
      }else if(!DE_by_cell){
        width=width+5
      }
      
      pdf(file=visFile, onefile = T, width=width, height=10)
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
          p0<-ggplot(geneexp, aes(x=Group, y=Gene, col=Group)) + geom_violin() + geom_jitter(width = 0.2) + facet_grid(~Sample) + theme_bw() + 
            scale_color_manual(values = groupColors) +
            NoLegend() + xlab("") + ylab("Gene Expression") + theme(strip.background=element_blank())
          
          p1<-DimPlot(cell_obj, reduction = "umap", label=T, group.by="DisplayGroup") + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5)) + xlim(xlim) + ylim(ylim)
          
          p2<-MyFeaturePlot(object = cell_obj, features=as.character(siggene), order=T)
          p<-p0+p1+p2+plot_layout(design="AA
BC")
          
          # p<-ggarrange(p0,                                                 # First row with scatter plot
          #              ggarrange(p1, p2, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
          #              nrow = 2, 
          #              labels = "A"                                        # Labels of the scatter plot
          # ) 
        }else{
          p0<-ggplot(geneexp, aes(x="1", y=Gene, color=Group)) + geom_violin() + geom_jitter(width = 0.2) + facet_grid(~Sample) + theme_bw() + 
            scale_color_manual(values = groupColors) +
            xlab("") + ylab("Gene Expression") + theme(strip.background=element_blank(), axis.text.x = element_blank())
          
          subcells<-colnames(cell_obj)[cell_obj$DisplayGroup == controlGroup]
          subobj<-subset(cell_obj, cells=subcells)
          p1<-MyFeaturePlot(object = subobj, features=siggene, order=T) + ggtitle(paste0("Control: ", controlGroup))
          p1<-suppressMessages(expr = p1 + xlim(xlim) + ylim(ylim) + fix.sc)
          
          subcells<-colnames(cell_obj)[cell_obj$DisplayGroup == sampleGroup]
          subobj<-subset(cell_obj, cells=subcells)
          p2<-MyFeaturePlot(object = subobj, features=siggene, order=T) + ggtitle(paste0("Sample: ", sampleGroup))
          p2<-suppressMessages(expr = p2  + xlim(xlim) + ylim(ylim) + fix.sc)
          
          if(!DE_by_cell){
            gd<-melt_cpm[melt_cpm$Gene==siggene,]
            g0<-ggplot(gd, aes(x=Group, y=CPM, color=Group)) + geom_violin() + geom_boxplot(width=0.2) + geom_jitter(width = 0.1) +  
              theme_bw3() +
              scale_color_manual(values = groupColors) +
              xlab("") + ylab("CPM") + NoLegend()
            p<-p0+g0+p1+p2+plot_layout(design="AAAAA
BCCDD")
          }else{
            p<-p0+p1+p2+plot_layout(design="AA
BC")
          }
        }
        p<-p+ plot_annotation(title=title)
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
  allcoords<-data.frame(obj@reductions$umap@cell.embeddings)
  allcoords$Cluster=obj[[cluster_name]]
  
  comp=unique(result$comparison)[1]
  for (comp in unique(result$comparison)) {
    compRes = result[result$comparison == comp,]
    rateMap=unlist(split(compRes$sigRate, compRes$cluster))
    
    obj$sigRate=rateMap[as.character(unlist(obj[[cluster_name]]))]
    
    pdf(paste0(outFile, ".", comp, ".sigGenePerc.pdf"), width=14, height=7)
    p1<-DimPlot(obj, group.by = cluster_name, label=T) + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5))
    p2<-MyFeaturePlot(obj, feature="sigRate", cols=c("lightgrey", "red")) + ggtitle("Percentage of DE genes in each cluster") + theme(plot.title = element_text(hjust=0.5))
    g<-ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))
    print(g)
    dev.off()
  }
}

write.csv(all_sigout, file=paste0(outFile, ".allsigout.csv"), quote=F)

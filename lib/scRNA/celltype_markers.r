source("scRNA_func.r")

library(Seurat)
library(ggplot2)
library(kableExtra)
library(dplyr)
library(scales)

if ("Seurat" %in% names(sessionInfo()$otherPkgs) & grepl("^4",sessionInfo()$otherPkgs$Seurat$Version)) { #Seurat version4
  logFcColName="avg_log2FC"
} else {
  logFcColName="avg_logFC"
}

options_table<-read.table("fileList1.txt", sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)

min.pct<-ifelse(is.null(myoptions$min.pct), 0.5, as.numeric(myoptions$min.pct))
logfc.threshold<-ifelse(is.null(myoptions$logfc.threshold), 0.6, as.numeric(myoptions$logfc.threshold))

finalList=readRDS(parFile1)
if(is.list(finalList)){
  all_obj=finalList$obj
  seurat_colors=finalList$seurat_colors
  celltype=read.csv(parFile3)
}else{
  all_obj=finalList
  celltype=all_obj@meta.data
  celltype<-unique(celltype[,c("seurat_clusters", myoptions$celltype_name)])
  colnames(celltype)<-c("seurat_clusters", "cell_type")
  seurat_colors<-hue_pal()(nrow(celltype))
}

draw_marker_genes<-function(all_obj, new.cluster.ids, file_prefix, celltype_prefix, min.pct, logfc.threshold ){
  assay=ifelse("SCT" %in% names(all_obj@assays), "SCT", "RNA")

  obj = subset(all_obj, seurat_clusters %in% names(new.cluster.ids))

  obj$cellactivity_clusters <- unlist(new.cluster.ids[as.character(unlist(obj$seurat_clusters))])
  
  clusterDf<-data.frame(seurat=unlist(obj[["seurat_clusters"]]), cellactivity=unlist(obj[["cellactivity_clusters"]]))
  clusterDf$seurat_colors<-seurat_colors[clusterDf$seurat]
  clusterDf$seurat_cellactivity<-paste0(clusterDf$seurat, " : ", clusterDf$cellactivity)
  seurat_cellactivity<-clusterDf$seurat_cellactivity
  
  clusterDf<-clusterDf[order(clusterDf$seurat),]
  caCount<-table(clusterDf$cellactivity)
  clusterDf$caCount<-caCount[clusterDf$cellactivity]
  clusterDf<-clusterDf[order(-clusterDf$caCount, clusterDf$seurat),]

  seurat_cellactivity<-factor(seurat_cellactivity, levels=unique(clusterDf$seurat_cellactivity))
  seurat_cellactivity_colors<-unique(clusterDf$seurat_colors)
  
  obj$seurat_cellactivity_clusters <-seurat_cellactivity
  
  clusters<-data.frame("cell" = c(1:length(obj$seurat_clusters)), "seurat_clusters"=as.numeric(as.character(obj$seurat_clusters)), "cellactivity_clusters"=obj$cellactivity_clusters, "seurat_cellactivity_clusters"=obj$seurat_cellactivity_clusters, stringsAsFactors = F)
  rownames(clusters)<-names(obj$seurat_clusters)
  write.csv(clusters, file=paste0(file_prefix, ".cluster.csv"))
  
  uniqueClusters=unique(clusters[,c("seurat_clusters", "cellactivity_clusters")])
  cellTypeClusters=tapply(uniqueClusters$seurat_clusters, uniqueClusters$cellactivity_clusters, list)
  
  ae=AverageExpression(obj, group.by = "seurat_clusters", assays = "SCT")$SCT
  aemax=colnames(ae)[max.col(ae, ties.method = "first")]
  genemax=data.frame(gene=rownames(ae), max_cluster=aemax)
  write.csv(genemax, paste0(file_prefix, ".genemax.csv"), row.names=F)
  
  allclusters=unique(clusters$seurat_clusters)
  
  all_bw_markers=NULL
  all_in_markers=NULL
  all_both_markers=NULL
  all_max_markers=NULL
  all_display_markers=NULL
  all_figures=NULL
  idx=1
  for(idx in c(1:length(cellTypeClusters))){
    ctc_name = names(cellTypeClusters)[idx]
    ctc_filename=gsub(" ", "_", ctc_name)
    cat(paste0("finding markers for ", ctc_name, "\n"))
    ctc=unlist(cellTypeClusters[idx])
    other=allclusters[!(allclusters %in% ctc)]
    c=2
    for(c in ctc){
      if(length(other) > 0){
        cat(paste0("  finding markers for cluster ", c, " between cell types\n"))
        bw_markers=find_markers(obj, by_sctransform=by_sctransform, ident.1=c, ident.2=other,logfc.threshold = logfc.threshold, min.pct = min.pct)
        cat(paste0("    ", nrow(bw_markers), " found between cell types\n"))
        
        bw_markers$cluster=c
        bw_markers$celltype=ctc_name
    
        all_bw_markers=rbind(all_bw_markers, bw_markers)
      }else{
        bw_markers=NA
      }
      suffix=""
      
      if(length(ctc) > 1){
        cat(paste0("  finding markers for cluster ", c, " in cell types\n"))
        other_ctc=ctc[ctc != c]
        in_markers=find_markers(obj, by_sctransform=by_sctransform, ident.1=c, ident.2=other_ctc,logfc.threshold = logfc.threshold, min.pct = min.pct)
        cat(paste0("    ", nrow(in_markers), " found in cell types\n"))
        in_markers$cluster=c
        in_markers$celltype=ctc_name
  
        all_in_markers=rbind(all_in_markers, in_markers)
        
        if(!is.na(bw_markers)){
          both_markers=bw_markers[rownames(bw_markers) %in% rownames(in_markers),,drop=F]
          cat(paste0("  there are ", nrow(both_markers), " markers found in both between and in cell types\n"))
          all_both_markers=rbind(all_both_markers, both_markers)
          
          if(nrow(both_markers) >= 5){
            cur_markers=both_markers
            suffix=".both"
          }else{
            cur_markers=bw_markers
            suffix=".between"
          }
        }else{
          cur_markers=in_markers
          suffix=".in"
        }
      }else{
        cur_markers=bw_markers
        suffix=".between"
      }
      c_filename=paste0(celltype_prefix, c, "_", ctc_filename, suffix, ".csv")
      write.csv(cur_markers, file=c_filename)
      
      cluster_max_genes=subset(genemax, genemax$max_cluster==as.character(c))
      max_markers=subset(cur_markers, rownames(cur_markers) %in% cluster_max_genes$gene)
      all_max_markers=rbind(all_max_markers, max_markers)
  
      if(nrow(max_markers) < 5){
        display_markers=cur_markers %>% top_n(n = 20, wt = .data[[logFcColName]])
      }else{
        display_markers=max_markers %>% top_n(n = 20, wt = .data[[logFcColName]])
        suffix=paste0(suffix, ".max")
      }
  
      cur_display_markers=rownames(display_markers)
      dot_filename=paste0(celltype_prefix, c, "_", ctc_filename, suffix, ".dot.pdf")
      pdf(file=dot_filename, width=14, height=7)
      g=DotPlot(obj, features=cur_display_markers, assay=assay, group.by="seurat_cellactivity_clusters" ) + 
        xlab("") + ylab("") + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust=1))
      print(g)
      dev.off()
      
      all_figures=rbind(all_figures, data.frame("Cluster"=c, "Dotfile"=dot_filename))
      all_display_markers=rbind(all_display_markers, data.frame("Cluster"=c, "Gene"=cur_display_markers))
    }
  }
  
  write.csv(all_display_markers, file=paste0(file_prefix, ".all_display_markers.csv"))
  write.csv(all_max_markers, file=paste0(file_prefix, ".all_max_markers.csv"))
  if(!is.null(all_bw_markers)){
    write.csv(all_bw_markers, file=paste0(file_prefix, ".all_bw_markers.csv"))
  }
  if(!is.null(all_in_markers)){
    write.csv(all_in_markers, file=paste0(file_prefix, ".all_in_markers.csv"))
  }
  write.csv(all_both_markers, file=paste0(file_prefix, ".all_both_markers.csv"))
  write.csv(all_figures, file=paste0(file_prefix, ".all_figures.csv"))
}

new.cluster.ids=split(celltype$cell_type, celltype$seurat_clusters)
draw_marker_genes(all_obj, new.cluster.ids, paste0(outFile, "_celltype"), "celltype_", min.pct=min.pct, logfc.threshold=logfc.threshold);

if("tcell_type" %in% colnames(celltype)){
  tcelltype=celltype[celltype$tcell_type != "",]
  new.cluster.ids=split(tcelltype$tcell_type, tcelltype$seurat_clusters)
  draw_marker_genes(all_obj, new.cluster.ids, paste0(outFile, "_tcelltype"), "tcelltype_");
}

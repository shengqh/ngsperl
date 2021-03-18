library(Seurat)
library(ggplot2)
library(kableExtra)
library(dplyr)

if ("Seurat" %in% names(sessionInfo()$otherPkgs) & grepl("^4",sessionInfo()$otherPkgs$Seurat$Version)) { #Seurat version4
  logFcColName="avg_log2FC"
} else {
  logFcColName="avg_logFC"
}

options_table<-read.table("fileList1.txt", sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)

finalList=readRDS(parFile1)
obj=finalList$obj

clusters=read_cell_cluster_file(parFile2)

uniqueClusters=unique(clusters[,c("seurat_clusters", "cellactivity_clusters")])
cellTypeClusters=tapply(uniqueClusters$seurat_clusters, uniqueClusters$cellactivity_clusters, list)

ae=AverageExpression(obj, group.by = "seurat_clusters", assays = "SCT")$SCT
aemax=colnames(ae)[max.col(ae, ties.method = "first")]
genemax=data.frame(gene=rownames(ae), max_cluster=aemax)
write.csv(genemax, paste0(outFile, ".genemax.csv"), row.names=F)

allclusters=unique(clusters$seurat_clusters)

all_bw_markers=NULL
all_in_markers=NULL
all_both_markers=NULL
all_max_markers=NULL
all_display_markers=NULL
all_figures=NULL
idx=12
for(idx in c(1:length(cellTypeClusters))){
  ctc_name = names(cellTypeClusters)[idx]
  ctc_filename=gsub(" ", "_", ctc_name)
  cat(paste0("finding markers for ", ctc_name, "\n"))
  ctc=unlist(cellTypeClusters[idx])
  other=allclusters[!(allclusters %in% ctc)]
  c=2
  for(c in ctc){
    cat(paste0("  finding markers for cluster ", c, " between cell types\n"))
    bw_markers=find_markers(obj, by_sctransform=by_sctransform, ident.1=c, ident.2=other)
    cat(paste0("    ", nrow(bw_markers), " found between cell types\n"))
    
    bw_markers$cluster=c
    bw_markers$celltype=ctc_name

    all_bw_markers=rbind(all_bw_markers, bw_markers)

    suffix=""
    
    if(length(ctc) > 1){
      cat(paste0("  finding markers for cluster ", c, " in cell types\n"))
      other_ctc=ctc[ctc != c]
      in_markers=find_markers(obj, by_sctransform=by_sctransform, ident.1=c, ident.2=other_ctc)
      cat(paste0("    ", nrow(in_markers), " found in cell types\n"))
      in_markers$cluster=c
      in_markers$celltype=ctc_name

      all_in_markers=rbind(all_in_markers, in_markers)
      
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
      cur_markers=bw_markers
      suffix=".between"
    }
    c_filename=paste0(c, "_", ctc_filename, suffix, ".csv")
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
    dot_filename=paste0(c, "_", ctc_filename, suffix, ".dot.pdf")
    pdf(file=dot_filename, width=14, height=7)
    g=DotPlot(obj, features=cur_display_markers, assay="SCT", group.by="seurat_cellactivity_clusters" ) + 
      xlab("") + ylab("") + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust=1))
    print(g)
    dev.off()
    
    all_figures=rbind(all_figures, data.frame("Cluster"=c, "Dotfile"=dot_filename))
    all_display_markers=rbind(all_display_markers, data.frame("Cluster"=c, "Gene"=cur_display_markers))
  }
}

write.csv(all_display_markers, file=paste0(outFile, ".all_display_markers.csv"))
write.csv(all_max_markers, file=paste0(outFile, ".all_max_markers.csv"))
write.csv(all_bw_markers, file=paste0(outFile, ".all_bw_markers.csv"))
write.csv(all_in_markers, file=paste0(outFile, ".all_in_markers.csv"))
write.csv(all_both_markers, file=paste0(outFile, ".all_both_markers.csv"))
write.csv(all_figures, file=paste0(outFile, ".all_figures.csv"))

# unique_cell_types=unique(clusters$cellactivity_clusters)
# 
# 
# clusterMarkers<-finalList$markers %>% group_by(cluster)
# write.csv(clusterMarkers, file=paste0(prefix, ".allmarkers.csv"), row.names=F, quote = F)
# 
# if ("Seurat" %in% names(sessionInfo()$otherPkgs) & grepl("^4",sessionInfo()$otherPkgs$Seurat$Version)) { #Seurat version4
#   logFcColName="avg_log2FC"
# } else {
#   logFcColName="avg_logFC"
# }
# 
# #top10 <- finalList$markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# top10 <- finalList$markers %>% group_by(cluster) %>% top_n(n = 10, wt = .data[[logFcColName]])
# top10marker_file = paste0(prefix, ".top10markers.csv")
# write.csv(top10, file=top10marker_file, row.names=F, quote = F)
# 
# top5 <- finalList$markers %>% group_by(cluster) %>% top_n(n = 5, wt = .data[[logFcColName]])
# print(kable(top5, caption=tabRef("top5", "Top 5 marker genes in each cluster")) %>% kable_styling() %>% htmltools::HTML())
# 
# fig.height=max(10, length(top10$gene) / 10)
# fig.width=12
# DoHeatmap(obj, features = top10$gene, group.colors=seurat_colors, angle = 90) + NoLegend()
# 
# fig.width=18
# fig.height=max(10, length(top10$gene) / 10)
# DoHeatmap(obj, features = top10$gene, group.by="seurat_cellactivity_clusters", group.colors=seurat_cellactivity_colors, angle=90) + NoLegend()
# ```


source("scRNA_func.r")

library(data.table)
library(Seurat)
library(heatmap3)

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

species=myoptions$species
markerfile<-myoptions$db_markers_file
remove_subtype<-myoptions$remove_subtype
annotate_tcell<-ifelse(myoptions$annotate_tcell == "0", FALSE, TRUE)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
assay=ifelse(myoptions$by_sctransform == "0", "RNA", "SCT")

remove_subtype_of=remove_subtype
if(annotate_tcell){
  if(remove_subtype_of != ""){
    remove_subtype_of = paste0(remove_subtype_of, ",T cells")
  }else{
    remove_subtype_of = "T cells"
  }
}

if(!file.exists(parFile1)){
  parFile1_tmp = gsub(".normByUpQuantile.csv", ".normByTotal.csv", parFile1)
  if(file.exists(parFile1_tmp)){
    parFile1 = parFile1_tmp
  }
}
data.norm=read.csv(parFile1, row.names=1,check.names = F)

cell_activity_database<-read_cell_markers_file(markerfile, species, remove_subtype_of, HLA_panglao5_file)
if("curated_markers_file" %in% names(myoptions)){
  curated_markerfile<-myoptions$curated_markers_file
  if (curated_markerfile != "") {
    curated_markers_df<-read.table(curated_markerfile, sep="\t", header=F, stringsAsFactors=F)
    curated_markers_celltype<-split(curated_markers_df$V2, curated_markers_df$V1)
    cellType=cell_activity_database$cellType
    for(cmct in names(curated_markers_celltype)){
      cellType[[cmct]]=curated_markers_celltype[[cmct]]
    }
    weight=calc_weight(cellType)
    cell_activity_database=list(cellType=cellType, weight=weight)
  }
}

predict_celltype<-ORA_celltype(data.norm,cell_activity_database$cellType,cell_activity_database$weight)

cta_ora_mat = get_cta_ora_mat(predict_celltype)
cell_activity_database$predicted<-predict_celltype
cell_activity_database$cta_mat=cta_ora_mat$cta_mat
cell_activity_database$ora_mat=cta_ora_mat$ora_mat

cell_type=list(cell_activity_database=cell_activity_database)

new.cluster.ids<-names(predict_celltype$max_cta)
names(new.cluster.ids) <- colnames(data.norm)

id_tbl=data.frame("seurat_clusters"=names(new.cluster.ids), "cell_type"=new.cluster.ids)
rownames(id_tbl)=id_tbl$seurat_clusters

if (annotate_tcell){
  tcell_activity_database<-get_selfdefinedCelltype(tcell_markers_file)
  tcell_clusters<-names(new.cluster.ids)[new.cluster.ids=="T cells"]
  tcell_data.norm<-data.norm[,tcell_clusters,drop=F]
  tcell_predict_celltype<-ORA_celltype(tcell_data.norm,tcell_activity_database$cellType,tcell_activity_database$weight)
  
  cta_ora_mat = get_cta_ora_mat(tcell_predict_celltype)
  tcell_activity_database$predicted<-tcell_predict_celltype
  tcell_activity_database$cta_mat=cta_ora_mat$cta_mat
  tcell_activity_database$ora_mat=cta_ora_mat$ora_mat
  
  tcell_new.cluster.ids<-names(tcell_predict_celltype$max_cta)
  names(tcell_new.cluster.ids) <- tcell_clusters
  
  id_tbl$tcell_type=""
  id_tbl[tcell_clusters,"tcell_type"]=tcell_new.cluster.ids
  
  new.cluster.ids[tcell_clusters] = tcell_new.cluster.ids
  
  cell_type$tcell_activity_database=tcell_activity_database
}

id_tbl$cellactivity_clusters=new.cluster.ids

cellColor <- function(style) 
{
  fg  <- style$getFillForegroundXSSFColor()
  rgb <- tryCatch(fg$getRgb(), error = function(e) NULL)
  rgb <- paste0("#", paste(rgb, collapse = ""))
  return(rgb)
}

if(file.exists(myoptions$summary_layer_file)){
  require(xlsx)
  layers=read.xlsx(myoptions$summary_layer_file, 1, header = TRUE)
  isna<-apply(layers, 2, function(x){
    all(is.na(x))
  })
  layers=layers[,!isna]
  
  wb     <- loadWorkbook(myoptions$summary_layer_file)
  sheet1 <- getSheets(wb)[[1]]
  
  # get all rows
  color_map = list()
  rows  <- getRows(sheet1)
  cells <- getCells(rows, colIndex = 6)
  for(i in c(2:length(cells))){
    v = getCellValue(cells[[i]])
    s = getCellStyle(cells[[i]])
    c = cellColor(s)
    color_map[[v]]=c
  }
  
  lastLayer=colnames(layers)[ncol(layers)]
  layers_map<-split(layers[, lastLayer], layers[,1])
  
  miss_celltype=id_tbl$cell_type[!(id_tbl$cell_type %in% names(layers_map))]
  for (mct in miss_celltype){
    layers_map[mct]=mct
  }
  id_tbl$summary_layer=unlist(layers_map[id_tbl$cell_type])
  
  saveRDS(color_map, file=paste0(outFile, ".summary_layer_color.rds"))
}

id_tbl$seurat_cellactivity_clusters = paste0(id_tbl$seurat_clusters, " : ", id_tbl$cell_type )
id_tbl$seurat_cellactivity_clusters=factor(id_tbl$seurat_cellactivity_clusters, levels=id_tbl$seurat_cellactivity_clusters)

write.csv(id_tbl, file=paste0(outFile, ".celltype.csv"), row.names=F)
saveRDS(cell_type, file=paste0(outFile, ".celltype.rds"))

png(paste0(outFile, ".score.png"), width=3000, height=3000, res=300)
heatmap3(cell_type$cell_activity_database$cta_mat, scale="none", margins=c(10,5))
dev.off()

clusters=read.csv(parFile2, header=T, stringsAsFactors = F, row.names=1)
clusters$cellactivity_clusters=new.cluster.ids[as.character(clusters$seurat_clusters)]
clusters$seurat_cellactivity_clusters=paste0(clusters$seurat_clusters, " : ", clusters$cellactivity_clusters)
write.csv(clusters, file=paste0(outFile, ".celltype_cluster.csv"))

draw_celltype_bar<-function(obj, seurat_clusters, celltypes, celltype_name, outFilePrefix){
  cmap<-split(seurat_clusters, celltypes)
  cts_clusters<-lapply(cmap, function(x){paste0(x, collapse = ",")})
  cts<-paste0(names(cmap), " (",cts_clusters ,")")
  names(cts)<-names(cmap)
  
  celltype_cluster<-unlist(cts[celltypes])
  ctmap<-split(celltype_cluster, id_tbl$seurat_clusters)
  
  meta_value<-unlist(ctmap[obj$seurat_clusters])
  
  tb<-table(meta_value)
  tb<-tb[order(tb)]
  meta_value<-factor(meta_value, levels=names(tb))
  
  tbl<-table(obj$orig.ident, meta_value)
  tbl<-tbl/rowSums(tbl)*100
  mtbl<-reshape2::melt(tbl)
  colnames(mtbl)<-c("Sample", "Celltype", "Perc")
  
  g<-ggplot(mtbl, aes(x=Sample, y=Perc, fill=Celltype)) + geom_bar(stat="identity") + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ylab("Percentage of total cell captured") + ggtitle(paste0("Frequency of ", celltype_name, "Cell Type per Sample"))
  
  png(paste0(outFilePrefix, ".bar.png"), width=2000, height=1400, res=300)
  print(g)
  dev.off()
}  

if(file.exists(parFile3)){
  library(ggplot2)
  finalList=readRDS(parFile3)
  obj=finalList$obj
  
  idmap = split(id_tbl$seurat_cellactivity_clusters, id_tbl$seurat_clusters)
  obj$seurat_cellactivity_clusters = unlist(idmap[as.character(obj$seurat_clusters)])
  
  cat("draw pictures ... ")
  
  draw_celltype_bar(obj, id_tbl$seurat_clusters, id_tbl$cell_type, "", paste0(outFile, ".celltype"))
  
  p1<-DimPlot(object = obj, reduction = 'umap', label=TRUE, group.by="seurat_cellactivity_clusters") + guides(colour = guide_legend(override.aes = list(size = 3), ncol=1)) + ggtitle("")
  p2<-DimPlot(object = obj, reduction = 'umap', label=FALSE, group.by="orig.ident") + ggtitle("")
  p=p1+p2
  png(paste0(outFile, ".cluster.png"), width=6600, height=3000, res=300)
  print(p)
  dev.off()
  
  png(paste0(outFile, ".celltype.label.png"), width=3300, height=3000, res=300)
  p1<-DimPlot(object = obj, reduction = 'umap', label=TRUE, group.by="seurat_cellactivity_clusters") + guides(colour = guide_legend(override.aes = list(size = 3), ncol=1)) + ggtitle("")
  print(p1)
  dev.off()
  
  png(paste0(outFile, ".celltype.nolabel.png"), width=3300, height=3000, res=300)
  p1<-DimPlot(object = obj, reduction = 'umap', label=FALSE, group.by="seurat_cellactivity_clusters") + guides(colour = guide_legend(override.aes = list(size = 3), ncol=1)) + ggtitle("")
  print(p1)
  dev.off()
  
  cat("Draw marker gene heatmap\n")
  obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
  write.csv(obj.markers, paste0(outFile, ".markers.csv"))
  obj@misc$markers<-obj.markers
  
  if('avg_log2FC' %in% colnames(obj.markers)){
    top10 <- obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  }else{
    top10 <- obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  }
  
  if (nrow(top10)>200) {
    genesize=5
  } else {
    if (nrow(top10)>100) {
      genesize=6
    }  else {
      genesize=7
    } 
  }
  
  Idents(obj)<-'seurat_cellactivity_clusters'
  png(paste0(outFile, ".heatmap.png"), width=6600, height=6000, res=300)
  print(DoHeatmap(obj, features = top10$gene,slot="data")+ theme(axis.text.y = element_text(size = genesize))) 
  dev.off()
  
  if(file.exists(myoptions$summary_layer_file)){
    cat("Draw figure in summary layers\n")
    idmap<-split(id_tbl$summary_layer, id_tbl$seurat_clusters)
    obj$summary_layer = unlist(idmap[as.character(obj$seurat_clusters)])
    
    draw_celltype_bar(obj, id_tbl$seurat_clusters, id_tbl$summary_layer, "Broad ", paste0(outFile, ".summary_layer"))
    
    summary_color=color_map[unique(obj$summary_layer)]
    
    p1<-DimPlot(object = obj, reduction = 'umap', label=TRUE, group.by="seurat_cellactivity_clusters") + 
      guides(colour = guide_legend(override.aes = list(size = 3), ncol=1)) +
      ggtitle("Cluster Cell Type")
    p2<-DimPlot(object = obj, reduction = 'umap', label=FALSE, group.by="summary_layer") + 
      scale_color_manual(values=summary_color) + 
      ggtitle("Broad Cell Type")
    p=p1+p2
    png(paste0(outFile, ".summary_layer.png"), width=6600, height=3000, res=300)
    print(p)
    dev.off()
    
    if(file.exists(parSampleFile2)){
      groups<-read.table(parSampleFile2, sep="\t", stringsAsFactors = F) 
      if(length(unique(groups$V2)) > 1){
        sample_map<-split(groups$V2, groups$V1)
        obj$group=unlist(sample_map[as.character(obj$orig.ident)])
        
        width=2000 * length(unique(groups$V2)) + 300
        png(paste0(outFile, ".summary_layer.group.png"), width=width, height=2000, res=300)
        p<-DimPlot(object = obj, reduction = 'umap', label=FALSE, group.by="summary_layer", split.by="group") + 
          ggtitle("") +
          scale_color_manual(values=summary_color)+
          annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
          annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
        print(p)
        dev.off()
        
        width=3000 * length(unique(groups$V2)) + 300
        png(paste0(outFile, ".celltype.group.label.png"), width=width, height=3000, res=300)
        p<-DimPlot(object = obj, reduction = 'umap', label=TRUE, group.by="seurat_cellactivity_clusters", split.by="group") + 
          guides(colour = guide_legend(override.aes = list(size = 3), ncol=1)) +
          ggtitle("") +
          annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
          annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
        print(p)
        dev.off()
        
        png(paste0(outFile, ".celltype.group.nolabel.png"), width=width, height=3000, res=300)
        p<-DimPlot(object = obj, reduction = 'umap', label=FALSE, group.by="seurat_cellactivity_clusters", split.by="group") + 
          guides(colour = guide_legend(override.aes = list(size = 3), ncol=1)) +
          ggtitle("") +
          annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
          annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
        print(p)
        dev.off()
      }
    }
  }
}

source("scRNA_func.r")

library(Seurat)
library(ggplot2)
library(ggpubr)
library(scales)

finalList<-readRDS(parFile1)

obj<-finalList$obj

groups_tbl<-read.table(parSampleFile2, sep="\t", stringsAsFactors = F)
groups=split(groups_tbl$V2, groups_tbl$V1)
obj$group = unlist(groups[obj$orig.ident])

ngroup=length(unique(groups_tbl$V2))

#using RNA assay for visualization
DefaultAssay(obj)<-"RNA"

filelist=NULL

genes_tbl<-read.table(parSampleFile1, sep="\t", stringsAsFactors = F)
gene=genes_tbl$V1[1]
for (gene in genes_tbl$V1){
  gene_list<-unlist(strsplit(gene, '/'))
  if(length(gene_list) != 2){
    stop(paste0(gene, " is not valid for gene ratio localization map."))
  }
  
  gene1 = gene_list[[1]]
  gene2 = gene_list[[2]]
  
  gdata<-FetchData(object = obj, gene_list)
  gdata$ratio = gdata[,gene1] - gdata[,gene2]
  coords<-data.frame(obj@reductions$umap@cell.embeddings)
  gdata<-cbind(coords, data.frame(group=obj$group), gdata)
  gdata1<-subset(gdata, gdata$ratio == 0)
  gdata2<-subset(gdata, gdata$ratio < 0)
  gdata3<-subset(gdata, gdata$ratio > 0)
  gdata_sort_1<-rbind(gdata1, gdata2, gdata3)
  gdata_sort_2<-rbind(gdata1, gdata3, gdata2)
  
  pngfile1 = paste0(outFile, ".", gene1, "_vs_", gene2, ".1.all.png")
  png(filename=pngfile1, width= 2400, height=2000, res=300)
  g<-ggplot(gdata_sort_1, aes(UMAP_1, UMAP_2, col=ratio)) + 
    geom_point(shape=1) +
    scale_color_gradient2(midpoint = 0, low = "green", mid = "gray",
                        high = "red", space = "Lab", limits=c(min(gdata$ratio), max(gdata$ratio))) +
    theme_bw() +
    theme(strip.background=element_rect(fill="white"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    ) + labs(color=paste0("log(", gene, ")"))
  print(g)
  dev.off()

  pngfile2 = paste0(outFile, ".", gene1, "_vs_", gene2, ".1.group.png")
  png(filename=pngfile2, width= ngroup * 2000 + 400, height=2000, res=300)
  g<-g + facet_grid(~group)
  print(g)
  dev.off()

  pngfile3 = paste0(outFile, ".", gene1, "_vs_", gene2, ".2.all.png")
  png(filename=pngfile3, width= 2400, height=2000, res=300)
  g<-ggplot(gdata_sort_2, aes(UMAP_1, UMAP_2, col=ratio)) + 
    geom_point(shape=1) +
    scale_color_gradient2(midpoint = 0, low = "green", mid = "gray",
                          high = "red", space = "Lab", limits=c(min(gdata$ratio), max(gdata$ratio))) +
    theme_bw() +
    theme(strip.background=element_rect(fill="white"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
    ) + labs(color=paste0("log(", gene, ")"))
  print(g)
  dev.off()
  
  pngfile4 = paste0(outFile, ".", gene1, "_vs_", gene2, ".2.group.png")
  png(filename=pngfile4, width= ngroup * 2000 + 400, height=2000, res=300)
  g<-g + facet_grid(~group)
  print(g)
  dev.off()

  filelist=rbind(filelist, data.frame("file"=paste0(getwd(), "/", c(pngfile1, pngfile2, pngfile3, pngfile4)), "genes"=paste0(gene1, "_vs_", gene2)))
}

write.csv(filelist, paste0(outFile, ".figure.files.csv"), row.names=F)
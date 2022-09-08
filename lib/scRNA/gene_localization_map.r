
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
  stop(paste0("There is missing genes:", paste0(miss_genes, ",")))
}
  
gene=genes_tbl$V1[1]
for (gene in genes_tbl$V1){
  gdata<-FetchData(object = obj, gene)
  colnames(gdata)<-"Gene"
  coords<-data.frame(obj@reductions$umap@cell.embeddings)
  gdata<-cbind(coords, data.frame(group=obj$group), gdata)
  gdata1<-subset(gdata, gdata$Gene == 0)
  gdata2<-subset(gdata, gdata$Gene > 0)

  g1<-ggplot(gdata, aes(UMAP_1, UMAP_2)) + 
    geom_point(data=gdata1, aes(col=Gene)) + 
    geom_point(data=gdata2, aes(col=Gene)) + 
    scale_colour_gradient(name=gene, low="grey", high="blue") + 
    theme_bw() +
    theme(strip.background=element_rect(fill="white"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
    )
  
  pngfile = paste0(outFile, ".", gene, ".all.png")
  png(filename=pngfile, width=1300, height=1100, res=300)
  print(g1)
  dev.off()
  filelist=rbind(filelist, data.frame("file"=paste0(getwd(), "/", pngfile), "gene"=gene, "type"="all"))

  if(has_group){
    g2<-g1+facet_grid(~group)
    pngfile = paste0(outFile, ".", gene, ".group.png")
    png(filename=pngfile, width= ngroup * 1200, height=1200, res=300)
    print(g2)
    dev.off()
    filelist=rbind(filelist, data.frame("file"=paste0(getwd(), "/", pngfile), "gene"=gene, "type"="group"))
  }
}

write.csv(filelist, paste0(outFile, ".figure.files.csv"), row.names=F)

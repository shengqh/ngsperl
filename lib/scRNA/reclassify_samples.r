library(reshape2)
library(ggplot2)
library(Seurat)
library(ggExtra)

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  h5file = "C:/Users/sheng/projects/paula_hurley/20201208_scRNA_split/filtered_feature_bc_matrix.h5"
  output_prefix="4701"
  cloupe_file="4701_VIP2_TvsN.csv"
  hashtag_regex = NA
}else{
  h5file = args[1]
  output_prefix = args[2]
  cloupe_file = args[3]
  hashtag_regex = NA
}

print(paste0("h5file=", h5file))
print(paste0("output_prefix=", output_prefix))
print(paste0("cloupe_file=", cloupe_file))
print(paste0("hashtag_regex=", hashtag_regex))

sdata<-Read10X_h5(h5file)

exp<-sdata[[1]]
meta<-sdata[[2]]
mat<-as.matrix(meta)
rowsum<-apply(mat>0, 1, sum)
mat<-mat[rowsum > (ncol(mat) / 2),]

if (!is.na(hashtag_regex)) {
  htos<-mat[grepl(hashtag_regex, rownames(mat)),]
}else{
  htos<-mat
}

tagnames=rownames(htos)

pbmc.hashtag <- CreateSeuratObject(counts = exp)
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

stypes=read.csv(cloupe_file, row.names=1)
pbmc.hashtag=subset(pbmc.hashtag, cells=rownames(stypes))
stypes=stypes[colnames(pbmc.hashtag),,drop=F]
pbmc.hashtag[["cloupe"]]=stypes[,1]

#object=pbmc.hashtag
#features=tagnames
#assay="HTO"
#identName="HTO.reclass"
rplot<-function(object, features, assay, identName, withAllCells=TRUE){
  DefaultAssay(object = object) <- assay
  data <- FetchData(object = object, vars = c(features, identName))
  mdata<-melt(data, id.vars=identName)
  if (withAllCells) {
    mdata2<-mdata
    mdata2[,1] = "All cells"
    mdata<-rbind(mdata, mdata2)
  }
  gfinal=list()
  for(feature in features){
    ddata=mdata[mdata$variable==feature,]
    g<-ggplot(ddata, aes_string(x="value")) + 
      geom_density() + 
      facet_grid(reformulate(".", identName), scale="free_y") + 
      xlab(feature) + theme_bw() + theme(strip.background=element_rect(colour="black", fill=NA))
    if (feature != features[1]){  
      g = g + ylab("")
    }
    gfinal = append(gfinal, list(g))
  }
  grid.arrange(grobs=gfinal, nrow=1)
}

png(paste0(output_prefix, "_cloupe.png"), width=2000, height=1600, res=300)
rplot(pbmc.hashtag, features, assay, "cloupe")
dev.off()

data <- FetchData(object=pbmc.hashtag, vars=c(tagnames, "cloupe"))
colnames(data)<-gsub("hto_","",colnames(data))

cloupe_names=unique(data$cloupe)
tagname=hto_tagnames[1]
cloupe_tag_map = NULL
for (tagname in tagnames) {
  tp=tapply(data[,tagname], data$cloupe, summary)
  tp_medians=lapply(tp, function(x){
    x["Median"]
  })
  tp_max=max(unlist(tp_medians))
  cloupe_name=names(tp_medians)[tp_medians==tp_max]
  
  negativeValues=data[,tagname][data$cloupe != cloupe_name]
  bp=boxplot(negativeValues,plot = FALSE)
  outlier=bp$stats[5]

  data[,cloupe_name] = ifelse(data[,tagname]>outlier, cloupe_name, "Negative")
  cloupe_tag_map = rbind(cloupe_tag_map, data.frame("cloupe"=cloupe_name, "tag"=tagname, "cutoff"=outlier))
}

data$Classification=apply(data, 1, function(x){
  xx=unique(x[cloupe_names])
  if (length(xx) > 1){
    xx = xx[xx != "Negative"]
    if (length(xx) > 1) {
      xx = "Douplet"
    }
  }
  return(xx)
})

tb=data.frame(table(data$Classification))
rownames(tb)=tb$Var1
tb$Rename=paste0(tb$Var1, " (", tb$Freq, ")")
data$ClassificationWithCount=tb[data$Classification,"Rename"]
pbmc.hashtag[["HTO.reclass"]]=data$Classification
pbmc.hashtag[["HTO.reclass_count"]]=data$ClassificationWithCount

png(paste0(output_prefix, "_reclass_count.png"), width=2000, height=1600, res=300)
rplot(pbmc.hashtag, features, assay, "HTO.reclass_count", withAllCells=FALSE)
dev.off()

if (length(hto_tagnames) == 2) {
  png(paste0(output_prefix, "_reclass_scatter.png"), width=2000, height=1600, res=300)
  p=ggplot(data, aes_string(x=tagnames[1], y=tagnames[2], color="ClassificationWithCount")) + geom_point() + 
    labs(color = "Sample") + theme_bw() + theme(legend.position="bottom")
  
  g<-ggMarginal(p,type = "histogram", margins = "both", size = 4, fill="white", xparams = list(bins=100), yparams=list(bins=100))  
  print(g)
  dev.off()
}

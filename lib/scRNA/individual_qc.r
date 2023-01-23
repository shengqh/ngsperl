library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RCurl)
library(knitr)
library(kableExtra)

source("scRNA_func.r")
source("markerCode_filter.R")

option_table<-read.table("fileList2.txt", sep="\t", stringsAsFactors = F)
myoptions<-split(option_table$V1, option_table$V2)
myoptions$mt_cutoff=as.numeric(myoptions$mt_cutoff)
myoptions$nCount_cutoff=as.numeric(myoptions$nCount_cutoff)
myoptions$nFeature_cutoff_max=as.numeric(myoptions$nFeature_cutoff_max)
myoptions$nFeature_cutoff_min=as.numeric(myoptions$nFeature_cutoff_min)
myoptions$pca_dims=as.numeric(myoptions$pca_dims)
myoptions$resolution=as.numeric(myoptions$resolution)
myoptions$Remove_MtRNA=is_one(myoptions$Remove_MtRNA)
myoptions$Remove_rRNA=is_one(myoptions$Remove_rRNA)
Mtpattern=myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
Remove_Mt_rRNA=myoptions$Remove_MtRNA
resolution=as.numeric(myoptions$resolution)

species=myoptions$species # Hs or Mm

#raw data file locations
SampleInfos<-read.table("fileList1.txt", stringsAsFactors = F)
colnames(SampleInfos)<-c("countfile", "SampleId")

if(file.exists('fileList3.txt')){
  hto_samples<-read.table('fileList3.txt', stringsAsFactors = F)
  hto_map=split(hto_samples$V1, hto_samples$V2)
  
  tag_tb<-read.table(myoptions$hto_sample_file, sep="\t", stringsAsFactors = F, header=T)
}else{
  hto_map=list()
  tag_tb=NULL
}

##cutoff dataframe for each sample to filter empty droplets
Cutoffs<-data.frame(nFeature_cutoff_min=myoptions$nFeature_cutoff_min ,
                    nFeature_cutoff_max=myoptions$nFeature_cutoff_max,
                    nCount_cutoff=myoptions$nCount_cutoff, 
                    mt_cutoff=myoptions$mt_cutoff, 
                    cluster_remove=c(""),stringsAsFactors = F)

#every sample share the same cutoff, otherwise put every parameter in the Cutoff dataframe
if (dim(Cutoffs)[1]==1){
  Cutoffs<-data.frame(lapply(Cutoffs, rep, nrow(SampleInfos)),stringsAsFactors = F)
}

transpose=FALSE

celltype_predictmethod="cta" # cta: cell activity; ora: over-represent

##specific marker database############

markerfile<-myoptions$markers_file

#these R will use markerfile, add_Subtype, species, so have to put at the end
add_Subtype<-FALSE

##compare to markerCode, add filter criteria for each sample

marker<-data.frame(fread(markerfile))
species_ind<-regexpr(species,marker[,1])
marker_species<-marker[species_ind>0 & marker$ubiquitousness.index<0.05,]
##change the gene symbol only keep the first letter capitalize
if (species=="Mm") {
  marker_species$official.gene.symbol<-paste0(substr(marker_species$official.gene.symbol,1,1),substr(tolower(marker_species$official.gene.symbol),2,nchar(marker_species$official.gene.symbol)))
}
cellType<-tapply(marker_species$official.gene.symbol,marker_species$cell.type,list)
##combine with other marker genes
if(file.exists(myoptions$curated_markers_file)){
  curated_markers_df<-read.table(myoptions$curated_markers_file, sep="\t", header=F, stringsAsFactors=F)
  curated_markers_celltype<-split(curated_markers_df$V2, curated_markers_df$V1)
  for(cmct in names(curated_markers_celltype)){
    cellType[[cmct]]=curated_markers_celltype[[cmct]]
  }
}

###if add crc-specific signatures # if add_subtype==NULL, only use subtype but not markerfile
if(is.null(add_Subtype)) {
  subtype<-read.csv(subtypefile,as.is=T,header=F)
  celltype_add<-tapply(subtype[,1],subtype[,2],list)
  cellType<-celltype_add
}else if (add_Subtype==TRUE) {
  subtype<-read.csv(subtypefile,as.is=T,header=F)
  celltype_add<-tapply(subtype[,1],subtype[,2],list)
  cellType<-c(cellType,celltype_add)
}

knownmarkers<- c("COL1A1","COL1A2","DCN","TAGLN","MYH11","ATCA2","PECAM1","CDH5","CLDN5","CD68","C1QB","LYZ2","NKAIN4","LGFBP5","CD248","WIF1","CHAD","SULT1D1","CD3D","CD3G","NKG7","CD209A","H2-AB1","H2-EB1","KCNA1","PLP1","MYL1","CD79A")

if(species=="Mm") {
  knownmarkers<-paste0(substr(knownmarkers,1,1),substr(tolower(knownmarkers),2,nchar(knownmarkers)))
}

object.list<-list()

i=1
for (i in 1:nrow(SampleInfos)) {
  SampleInfo<-SampleInfos[i,]
  Cutoff<-Cutoffs[i,]
  bubble_file=myoptions$bubblemap_file
  Ensemblfile=NULL
  info<-preprocess(SampleInfo,
                     Cutoff,
                     Mtpattern,
                     resolution, 
                     Remove_Mt_rRNA,
                     celltype_predictmethod,
                     transpose=transpose,
                     hto_map=hto_map,
                     tag_tb=tag_tb,
                     Ensemblfile=Ensemblfile,
                     bubble_file=bubble_file)
  
  object.list<-c(object.list, info)
}

save(object.list,file=paste0("objectlist.Rdata"))
#load("objectlist.Rdata")

stats<-lapply(object.list, function(x){unlist(x$preprocess)})
stats_df<-data.frame(do.call(rbind, stats))
colnames(stats_df)<-gsub("preprocess.","",colnames(stats_df))
write.table(stats_df, file="qc_filter_config.txt", sep="\t", row.names=F)

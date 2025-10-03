
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sf)
library(Matrix)
library(cowplot)
library(SpatialExperiment)


marker<-data.frame(fread(markerfile))
species_ind<-regexpr(species,marker[,1])
marker_species<-marker[species_ind>0 & marker$ubiquitousness.index<0.05,]
# bubble_file <- data.frame(read.delim(bubblemap_file,header=F))
##change the gene symbol only keep the first letter capitalize
if (species=="Mm") {
  marker_species$official.gene.symbol<-paste0(substr(marker_species$official.gene.symbol,1,1),substr(tolower(marker_species$official.gene.symbol),2,nchar(marker_species$official.gene.symbol)))
}
##
cellType<-tapply(marker_species$official.gene.symbol,marker_species$cell.type,list)


###if add crc-specific signatures # if add_subtype==NULL, only use subtype but not markerfile

if(is.null(add_Subtype)) {
  subtype<-read.csv(markerfile,as.is=T,header=F)
  #celltype_add<-tapply(subtype[,2],subtype[,1],list)
  #cellType<-celltype_add
} else if (add_Subtype==TRUE) {
  subtype<-read.csv(subtypefile,header=TRUE,sep="\t")
  #subtype<-subtype[subtype$source=="Wells"|subtype$source=="Bornstein",]
  celltype_add<-tapply(subtype[,3],subtype[,1],list)
  celltype_add2 <- lapply(celltype_add,function(x) {x[1:50]})
  cellType<-c(cellType,celltype_add2)
}

################
# Prepare Cell Cycle Scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

if (species=="Mm") {
  s.genes<-paste0(substr(s.genes,1,1),substr(tolower(s.genes),2,nchar(s.genes)))
  g2m.genes<-paste0(substr(g2m.genes,1,1),substr(tolower(g2m.genes),2,nchar(g2m.genes)))
}

################

ORA_celltype<-function(medianexp,cellType,method="cta"){
  freq<-sort((table(unlist(cellType)))/length(cellType))
  weight<-1+sqrt((max(freq)-freq)/(max(freq)-min(freq)))	
  
  ORA_result<-matrix(NA, nrow=length(cellType),ncol=dim(medianexp)[2])
  CTA_result<-matrix(0,nrow=length(cellType),ncol=dim(medianexp)[2])
  exp_z<-scale(medianexp)
  genenames<-rownames(medianexp)   
  for (j in 1: dim(medianexp)[2]){
    clusterexp<-medianexp[,j] 
    clusterexp_z<-exp_z[,j]
    for (i in 1:length(cellType)){
      ct_exp_names<-(intersect(genenames[clusterexp>0],cellType[[i]]))
      ct_exp<-length(ct_exp_names)
      ct_not_exp<-length(cellType[[i]])-ct_exp
      exp_not_ct<-sum(clusterexp>0)-ct_exp
      not_exp_not_ct<-length(clusterexp)-ct_not_exp 
      cont.table<-matrix(c(ct_exp,ct_not_exp,exp_not_ct,not_exp_not_ct),nrow=2)
      ORA_result[i,j]<-fisher.test(cont.table,alternative="greater")$p.value
      ###
      weight_ss<-weight[names(weight)%in%cellType[[i]]]
      ind<-match(names(weight_ss),genenames)
      exp_ss<-clusterexp_z[ind[!is.na(ind)]]
      weight_ss<-weight_ss[!is.na(ind)]
      CTA_result[i,j]<-sum(exp_ss*weight_ss)/(length(exp_ss)^0.3)
    }
  }
  rownames(ORA_result)<-rownames(CTA_result)<-names(cellType)
  minp_ora_ind<- apply(ORA_result,2,function(x){which.min(x)})
  minp_ora<-apply(ORA_result,2,min)
  names(minp_ora)<-rownames(ORA_result)[minp_ora_ind]
  
  max_cta_ind<- apply(CTA_result,2,function(x){which.max(x)})
  max_cta<-apply(CTA_result,2,max,na.rm=T)
  names(max_cta)<-rownames(CTA_result)[max_cta_ind]
  
  if (method=="cta") {
    cta_val<-max_cta
    p_val<-ORA_result[cbind(match(names(max_cta),rownames(ORA_result)),seq(1:length(max_cta)))]
    
    predict_result<-cbind(cta_val,p_val)
  }
  if (method=="ora") {
    p_val<-minp_ora
    cta_val<-CTA_result[cbind(match(names(minp_ora),rownames(CTA_result)),seq(1:length(minp_ora)))]
    
    predict_result<-cbind(cta_val,p_val)  
  }
  
  return(list(predict_result=predict_result, method=method,ora=ORA_result,cta=CTA_result,minp_ora=minp_ora,max_cta=max_cta,ct_exp_names=ct_exp_names))
}


PurpleAndYellow <- function(k = 50) {
  return(CustomPalette(low = "magenta", high = "yellow", mid = "black", k = k))
}

##plot the heatmap for marker genes in each cell type, annotate marker genes based on the existence of the database
Doheatmap_cellmarker_cellType<-function(v3d,top10,cellType,predict_celltype) {
  genemarker<-rev(top10$gene)
  orderind<-order(v3d@active.ident)
  data_m<-GetAssayData(v3d,slot="data")
  geneind<-match(genemarker,rownames(data_m))
  data_h<-data_m[geneind,orderind]
  data_h_scale<-scale(t(as.matrix(data_h)))
  data_h_scale<-MinMax(data_h_scale,min=-2.5,max=2.5)
  
  matchid<-NULL
  ###
  for (i in 1:length(levels(v3d@active.ident))){
    
    matchid<-c(matchid,ifelse(is.na(match(top10$gene[top10$cluster==i-1],cellType[[rownames(predict_celltype$predict_result)[i]]])),0,1))
    
  }
  matchid<-rev(matchid)
  predict_result<-data.frame(predict_celltype$predict_result)
  
  rowrange<-seq(from = 0, to = 1, length =ncol(data_h_scale))
  colrange<-seq(from = 0, to = 1, length =nrow(data_h_scale))
  clustnum<-table(v3d@active.ident)
  colind<-c(0,cumsum(clustnum))+c(clustnum/2,0)
  colind<-floor(colind)[-length(colind)]
  
  col_text<-ifelse(predict_result$p_val<1e-5,"blue","gray")
  names(colind)<-levels(v3d@active.ident)
  
  markernum<-table(tapply(top10$gene,top10$cluster))
  markerind<-cumsum(rev(markernum))
  
  old.par <- par(no.readonly = TRUE)  
  par(mar=c(5,5,8,1))
  image(data_h_scale,col=PurpleAndYellow(),axes=F)
  abline(v=colrange[cumsum(clustnum)],col="white",lwd=2)
  abline(h=rowrange[markerind]+0.5/length(rowrange),col="white",lwd=2)
  axis(side=2,at= rowrange[matchid==0],labels=genemarker[matchid==0],las=1,tick=FALSE,cex.axis=0.5,col.axis="red")
  axis(side=2,at= rowrange[matchid==1],labels=genemarker[matchid==1],las=1,tick=FALSE,cex.axis=0.5,col.axis="black")
  text(colrange[colind],1+40/length(colrange), labels=names(colind),srt=45,xpd=T,cex=0.8,col=col_text)
  text(colrange[colind],-20/length(colrange),labels=paste0("cta:",round(predict_result$cta_val,1)),srt=45,xpd=T,cex=0.8) 
  text(colrange[colind+15],-20/length(colrange),labels=paste0("p:",format(predict_result$p_val,format="e",digits=1)),srt=45,xpd=T,cex=0.8)
  par(old.par)
  
}


###plot the predict result 
Plot_predictcelltype<-function(predict_celltype,method="cta",topn=3) {
  cta_index<-apply(predict_celltype$cta,2,function(x){return(order(x,decreasing=T)[1:topn])})
  cta_index<-unique(sort(cta_index))
  cta_mat<- predict_celltype$cta[cta_index,]
  
  ora_index<-apply(predict_celltype$ora,2,function(x){return(order(x,decreasing=F)[1:topn])})
  ora_index<-unique(sort(ora_index))
  ora_mat<- predict_celltype$ora[ora_index,]
  ora_mat<--log10(ora_mat) 
  if (method=="cta") plotfunction_celltype(cta_mat)
  if (method=="ora") plotfunction_celltype(ora_mat)
  
  
}

##
plotfunction_celltype<-function(cta_mat) {
  
  old.par <- par(no.readonly = TRUE)  
  
  layout(matrix(c(2,0,0,1), 2, 2, byrow = TRUE),width=c(1,6),heights=c(1,6))
  par(mar=c(3,1,1,10))
  par(xaxs="i")
  par(yaxs="i")
  image(1:dim(cta_mat)[2],1:dim(cta_mat)[1],t(cta_mat),col=colorRampPalette(c("white","red"))(100),axes=F,xlab="",ylab="")
  axis(4,1:dim(cta_mat)[1],labels=rownames(cta_mat),las=2,tick=F,cex.axis=0.8)
  axis(1,1:dim(cta_mat)[2],labels=0:(dim(cta_mat)[2]-1),las=2,tick=F,cex.axis=0.8)
  
  par(mar=c(2,1,2,1))
  maxval<-round(max(cta_mat),digits=1)
  zval<-seq(0,maxval,0.1)
  image(1:length(zval),1,matrix(zval,ncol=1),col= colorRampPalette(c("white", "red"))(100),axes=F,xlab="",ylab="")
  mtext("0",side=1,line=0.5,at=1,cex=0.8)
  mtext(maxval,side=1,line=0.5,at=length(zval),cex=0.8)
  par(old.par)
  
}  


## if provide ensembl ID, using the convertfile (csv format) to convert to gene symbol, the first line is gene symbol and then second is ensembl
convertEnsembltoGeneSymbol<-function(counts,convertfile) {
  
  converttable<-read.csv(convertfile,stringsAsFactors = F)
  ind<-match(rownames(counts),converttable[,2])
  if (sum(is.na(ind)>0)) {
    counts<-counts[!is.na(ind),]
    genename<-converttable[ind[!is.na(ind)],1]
    warning('convertfile is incorrect')} else {genename<-converttable[ind,1]}
  
  
  ind_dup<-duplicated(genename)
  counts<-counts[!ind_dup,]
  rownames(counts)<-genename[!ind_dup]
  return(counts)
}



savelevel3and4<-function(SCLC,SampleInfo,Cutoff,resolution,outputdir) {
  
  SampleId<-SampleInfo$SampleId
  if (!dir.exists(paste0(outputdir, SampleId))) {
    if (!dir.create(paste0(outputdir, SampleId)))
    {stop("folder cannot be created")}
  }
  outputdir<-paste0(outputdir,SampleId)
  
  write.csv(c(SampleInfo,Cutoff,UMI_median=median(SCLC@meta.data$nCount_RNA),nGene_median=median(SCLC@meta.data$nFeature_RNA),nCells=ncol(SCLC),resolution=resolution),file=paste0(outputdir,"/QCparam.csv"))
  
  countdata<-GetAssayData(SCLC,slot="counts")
  normalizeddata<-  GetAssayData(SCLC,slot="data")
  pca<-Embeddings(SCLC,reduction="pca")[,1:2]
  umap<-Embeddings(SCLC,reduction="umap")
  
  write.csv(countdata,file=paste0(outputdir,"/filteredcounts.csv"))
  write.csv(normalizeddata,file=paste0(outputdir,"/normalizedcounts.csv"))
  write.csv(data.frame(pca,umap,Cell_ann=SCLC@active.ident),file=paste0(outputdir,"/cellannotation.csv"))
  
}

getMedianExp<-function(SCLC, chunk=1000) {
  
  data<-GetAssayData(SCLC,slot="data",assay="RNA")
  ngenes <- nrow(data)
  ncol<-ncol(data)
  
  if (ngenes > chunk) {
    by.chunk <- cut(seq_len(ngenes), ceiling(ngenes/chunk))
  } else {
    by.chunk <- factor(integer(ngenes))
  }
  
  medianExp<-matrix(NA,nrow=ngenes,ncol=length(unique(SCLC@active.ident)))
  for (element in levels(by.chunk)) {
    current <- by.chunk == element
    cur.exprs <- data[current, , drop = FALSE]
    medianExp[current,]<-t(apply(cur.exprs,1,function(x){tapply(x,SCLC@active.ident,median)}))
    
  }
  rownames(medianExp)<-rownames(data)
  
  
  return(medianExp)
}

getMeanExp<-function(SCLC, chunk=1000) {
  
  data<-GetAssayData(SCLC,slot="data",assay="RNA")
  ngenes <- nrow(data)
  ncol<-ncol(data)
  
  if (ngenes > chunk) {
    by.chunk <- cut(seq_len(ngenes), ceiling(ngenes/chunk))
  } else {
    by.chunk <- factor(integer(ngenes))
  }
  
  meanExp<-matrix(NA,nrow=ngenes,ncol=length(unique(SCLC@active.ident)))
  for (element in levels(by.chunk)) {
    current <- by.chunk == element
    cur.exprs <- data[current, , drop = FALSE]
    meanExp[current,]<-t(apply(cur.exprs,1,function(x){tapply(x,SCLC@active.ident,mean)}))
    
  }
  rownames(meanExp)<-rownames(data)
  
  
  return(meanExp)
}

getSumExp<-function(SCLC, chunk=1000) {
  
  data<-GetAssayData(SCLC,slot="counts")
  ngenes <- nrow(data)
  
  
  if (ngenes > chunk) {
    by.chunk <- cut(seq_len(ngenes), ceiling(ngenes/chunk))
  } else {
    by.chunk <- factor(integer(ngenes))
  }
  
  SumExp<-NULL
  for (element in levels(by.chunk)) {
    current <- by.chunk == element
    cur.exprs <- data[current, , drop = FALSE]
    groups<-data.frame(SCLC[["PatientID"]], SCLC[["Location"]],SCLC[["Histo"]], SCLC[["Disease"]],SCLC[["cellcluster_ann"]],stringsAsFactors=F) 
    tmp<-cbind(groups,t(as.matrix(cur.exprs)))
    tmp_exp<- tmp %>% group_by(PatientID, Location, Histo, Disease,cellcluster_ann) %>% summarise_each(funs(sum))
    SumExp<-cbind(SumExp,as.matrix(tmp_exp[,-c(1:5)]))
    if (element==levels(by.chunk)[1]) {
      groupInfo<-as.data.frame(tmp_exp[,c(1:5)]) 
    }
  }
  return(list(Exp=t(SumExp),group=groupInfo))
}


get_bubble_plot<-function(obj, 
                          cur_res, 
                          cur_celltype, 
                          bubblemap_file, 
                          assay="RNA", 
                          orderby_cluster=FALSE, 
                          split.by=NULL, 
                          rotate.title=TRUE, 
                          group.by=NULL, 
                          use_blue_yellow_red=TRUE, 
                          species="Hs",
                          dot.scale=6){
  
  old_assay = DefaultAssay(obj)
  DefaultAssay(obj) = assay
  allgenes=rownames(obj)
  DefaultAssay(obj) = old_assay
  
  if(assay=="ADT"){
    genes_df <- read_bubble_genes(bubblemap_file, allgenes, species=NULL)
    if(nrow(genes_df) == 0){
      stop(paste0("no proteins in ", bubblemap_file," found in data, please double check your protein names. They should match with following protein names: ", paste0(allgenes, collapse = ", ")))
    }
  }else{
    genes_df <- read_bubble_genes(bubblemap_file, allgenes, species=species)
    if(nrow(genes_df) == 0){
      stop(paste0("no genes in ", bubblemap_file," found in data, is it posible wrong species used in function? species=", species))
    }
  }
  
  gene_groups=split(genes_df$gene, genes_df$cell_type)
  
  if(is.null(group.by)){
    if(is.null(cur_celltype)){
      stop("cur_celltype cannot be null in get_bubble_plot if group.by is null.")
    }
    if(is.null(cur_res) || is.na(cur_res)){
      cur_res = paste0(cur_celltype, "_cluster")
      obj<-build_dummy_cluster(obj, label.by=cur_celltype, new_cluster_name=cur_res)
    }
    
    cell_type=obj@meta.data
    cell_type$cell_type <- cell_type[,cur_celltype]
    
    # ct_levels<-c("B cells", "Plasma cells", "NK cells", "T cells", "Macrophages", "Dendritic cells", "Monocytes", "Mast cells", "Endothelial cells", "Fibroblasts", "Epithelial cells", "Basal cells", "Olfactory epithelial cells", "Ciliated cells")
    # ct<-cell_type[!duplicated(cell_type$cell_type),]
    # missed = ct$cell_type[!(ct$cell_type %in% ct_levels)]
    # if(length(missed) > 0){
    #   ct_levels = c(ct_levels, as.character(missed))
    # }
    # ct_levels = ct_levels[ct_levels %in% ct$cell_type]
    # cell_type$cell_type<-factor(cell_type$cell_type, levels=ct_levels)
    
    if (! is.null(split.by)){ 
      if(orderby_cluster){
        cell_type<-cell_type[order(cell_type[,cur_res], cell_type[,split.by]),]
      }else{
        cell_type<-cell_type[order(cell_type$cell_type, cell_type[,cur_res], cell_type[,split.by]),]
      }
      cell_type$seurat_celltype_clusters=paste0(cell_type[,cur_res], ": ", cell_type$cell_type, ": ", cell_type[,split.by])
    }else{
      if(orderby_cluster){
        cell_type<-cell_type[order(cell_type[,cur_res]),]
      }else{
        cell_type<-cell_type[order(cell_type$cell_type, cell_type[,cur_res]),]
      }
      cell_type$seurat_celltype_clusters=paste0(cell_type[,cur_res], ": ", cell_type$cell_type)
    }
    cell_type$seurat_celltype_clusters=factor(cell_type$seurat_celltype_clusters, levels=unique(cell_type$seurat_celltype_clusters))
    group.by="seurat_celltype_clusters"
    
    cell_type<-cell_type[colnames(obj),]
    obj@meta.data<-cell_type
  }
  
  g<-get_dot_plot(obj, group.by, gene_groups, assay, rotate.title=rotate.title, use_blue_yellow_red=use_blue_yellow_red, dot.scale=dot.scale)
  
  return(g)
}


get_dot_plot<-function(obj, group.by, gene_groups, assay="RNA", rotate.title=TRUE, use_blue_yellow_red=TRUE, dot.scale=6){
  genes=unique(unlist(gene_groups))
  assaydata=MyGetAssayData(obj, assay=assay, slot="data")
  if(!all(genes %in% rownames(assaydata))){
    missed_genes = genes[!(genes %in% rownames(assaydata))]
    missed_genes=missed_genes[c(1:min(5, length(missed_genes)))]
    stop(paste0("some genes are not in ", assay, " assay, here is the first few:", paste0(missed_genes, collapse = ",")))
  }
  
  g<-DotPlot(obj, features=genes, assay=assay, group.by=group.by)
  
  gdata<-g$data
  
  data.plot<-NULL
  gn=names(gene_groups)[1]
  for(gn in names(gene_groups)){
    gs=gene_groups[[gn]]
    gdd<-gdata[gdata$features.plot %in% gs,]
    if(nrow(gdd)== 0){
      stop(gn)
    }
    gdd$feature.groups=gn
    data.plot<-rbind(data.plot, gdd)
  }
  
  data.plot$feature.groups=factor(data.plot$feature.groups, levels=names(gene_groups))
  
  color.by <- "avg.exp.scaled"
  scale.func <- scale_radius
  scale.min = NA
  scale.max = NA
  cols = c("lightgrey", "blue")
  
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "id")) + 
    geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    guides(size = guide_legend(title = "Percent Expressed")) +
    labs(x = "Features", y = "Identity") +
    theme_cowplot() + 
    # facet_grid(facets = ~feature.groups, scales = "free_x", space = "free_x", switch = "y") + 
    facet_grid(rows = ~feature.groups, scales = "free_x", space = "free_x", switch = "y") + 
    theme(panel.spacing = unit(x = 1,units = "lines"), strip.background = element_blank())
  
  if(use_blue_yellow_red){
    plot <- plot + scale_colour_gradient2(low="blue", mid="yellow", high="red", midpoint=0 )
  }else{
    plot <- plot + scale_color_gradient(low="lightgray", high="blue")
  }
  
  g=plot + xlab("") + ylab("") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                    axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
                                                    strip.background = element_blank())
  
  if(rotate.title){
    g=g+theme(strip.text.x = element_text(angle=90, hjust=0, vjust=0.5))
  }
  return(g)
}

read_bubble_genes<-function(bubblemap_files, allgenes=c(), species="Hs"){
  result = NULL
  bubblemap_file=bubblemap_files[1]
  for(bubblemap_file in bubblemap_files){
    genes = do_read_bubble_genes(bubblemap_file, allgenes, species)
    if(is.null(result)){
      result = genes
    }else{
      result = rbind(result, genes)
    }
  }
  result$cell_type=factor(result$cell_type, levels=unique(result$cell_type))
  
  return(result)
}

do_read_bubble_genes<-function(bubblemap_file, allgenes=c(), species="Hs"){
  library("readxl")
  library("tidyr")
  
  if(grepl(".txt$", bubblemap_file)){
    genes <- read.table(bubblemap_file, header=F, sep="\t", stringsAsFactors = F)
    if(any(grepl(",", genes$V2))){
      genes <- genes[,c(2,1)]
    }
  }else{
    genes <- data.frame(read_xlsx(bubblemap_file, sheet = 1))
    if(colnames(genes)[1] %in% c('Count', 'Index')){
      genes <- genes[,c(2:ncol(genes))]
    }
    
    while(all(is.na(genes[,1]))){
      genes <- genes[,c(2:ncol(genes))]
    }
    
    if(colnames(genes)[2] == 'Marker.Gene'){
      genes <- genes[,c(2:ncol(genes))]
    }
  }
  
  colnames(genes)[1:2] = c("gene", "cell_type")  
  
  #don't use the marker genes other than thhe first column, it might just annotation
  #if you want to use those genes, put into the first column
  # if(ncol(genes) > 2){
  #   #for some excel files, the third column is addtional genes
  #   if(any(grepl(',', genes[,1])) & any(grepl(',', genes[,3]))){
  #     genes[,3][is.na(genes[,3])]<-""
  #     genes[,1] = paste0(genes[,1], ",", genes[,3])
  #   }
  
  #   genes<-genes[,c(1,2)]
  # }
  
  for(idx in c(2:nrow(genes))){
    if(is.na(genes[idx,"cell_type"])){
      genes[idx,"cell_type"]=genes[idx-1,"cell_type"]
    }
  }
  
  #for some excel file, multiple genes are in same row but seperated by ',', 
  #use sepearte_rows to put them in different rows.
  genes<-data.frame(separate_rows(genes, gene, sep="[, ]+"))
  genes<-genes[genes$gene != "",]
  
  if(!is.null(species)){
    if(tolower(species) == "mm" | tolower(species) == "mouse"){
      genes$gene = toMouseGeneSymbol(genes$gene)
    }
  }
  
  gene_names=genes$gene
  gene_names[gene_names=="PECAM"] = "PECAM1"
  gene_names[gene_names=="HGD1B"] = "HGD"
  gene_names[gene_names=="EpCAM"] = "EPCAM"
  gene_names[gene_names=="CD25"] = "IL2RA"
  gene_names[gene_names=="ACTAA2"] = "ACTA2"
  gene_names[gene_names=="MTND6"] = "MT-ND6"
  gene_names[gene_names=="FOXJ!"] = "FOXJ1"
  
  genes$gene<-gene_names
  
  genes<-genes[!duplicated(genes),]
  
  if(length(allgenes) > 0){
    miss_genes=setdiff(genes$gene, allgenes)
    writeLines(miss_genes, con="miss_gene.csv")
    
    genes<-genes[genes$gene %in% allgenes,]
  }
  genes$cell_type=factor(genes$cell_type, levels=unique(genes$cell_type))
  
  return(genes)
}

build_dummy_cluster<-function(obj, label.by, new_cluster_name, new_cluster_name_label=paste0(new_cluster_name, "_label")){
  groups<-as.character(obj@meta.data[,label.by])
  gt<-table(groups)
  gt<-gt[order(gt, decreasing=T)]
  dummy_cluster<-c(0:(length(gt)-1))
  names(dummy_cluster)<-names(gt)
  dc<-factor(dummy_cluster[groups], levels=dummy_cluster)
  obj@meta.data[,new_cluster_name]<-dc
  obj@meta.data[,new_cluster_name_label]<-paste0(obj@meta.data[,new_cluster_name], ": ", groups)
  return(obj)
}

MyGetAssayData<-function(obj, assay, slot){
  if(is_assay_5_plus(obj, assay)){
    cur_assay=obj[[assay]]
    return(LayerData(cur_assay, layer=slot))
  }else{
    return(GetAssayData(obj, assay=assay, slot=slot))
  }
}

is_assay_5_plus<-function(obj, assay){
  return(class(obj[[assay]]) == 'Assay5')
}


seurat_to_spe <- function(seu, sample_id, img_id) {
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## Extract spatial coordinates
  spatialCoords <- as.matrix(
    seu@images[[img_id]]@boundaries$centroids@coords[, 1:2])
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  spe
}

fill_symmetric_na <- function(mat) {
  is_na <- is.na(mat)
  mat[is_na] <- t(mat)[is_na]
  return(mat)
}

plot_type1_transformation <- function(fit_dat_sub, stats) {
  plot_dat <- fit_dat_sub[, c("celltype1", "celltype2", stats)]
  cts <- unique(c(plot_dat$celltype1, plot_dat$celltype2))
  sup_dat <- data.frame(
    "celltype1" = cts,
    "celltype2" = cts
  )
  sup_dat[, stats] <- 0
  plot_dat <- rbind(plot_dat, sup_dat)
  plot_dat <- reshape(plot_dat,
                      idvar = "celltype1",
                      timevar = "celltype2",
                      direction = "wide"
  ) |>
    as.data.frame(optional = TRUE) |>
    col2rownames("celltype1")
  colnames(plot_dat) <- gsub(
    paste0("^", stats, "\\."), "",
    colnames(plot_dat)
  )
  plot_dat <- plot_dat[colnames(plot_dat), ]
  filled_data <- fill_symmetric_na(plot_dat)
  return(filled_data)
}

plot_type2_transformation <- function(fit_dat, stats) {
  fit_dat <- as.data.frame(fit_dat[, c(
    "celltype1", "celltype2",
    "ROI", stats
  )])
  fit_dat[["pair"]] <- paste0(
    fit_dat[["celltype1"]], "-",
    fit_dat[["celltype2"]]
  )
  fit_dat <- fit_dat[, c("pair", "ROI", stats)]
  filled_data <- reshape(fit_dat,
                         idvar = "pair", timevar = "ROI",
                         direction = "wide"
  ) |> col2rownames(column=1)
    # col2rownames("pair")
  colnames(filled_data) <- gsub(
    paste0("^", stats, "\\."), "ROI#",
    colnames(filled_data)
  )
  return(filled_data)
}

plotCorHeatmap_modified <- function (model.result, stats = c("cor.coef", "t", "p.Pos", "p.Neg"), 
          roi = "all", cell.type = "all") 
{
  if (!all(c("cor.coef", "p.Pos", "p.Neg") %in% colnames(model.result))) {
    stop("Please run corDensity before using this function.")
  }
  fit_dat <- model.result
  if (length(stats) != 1) {
    stats <- "cor.coef"
  }
  else if (!(stats %in% c("cor.coef", "t", "p.Pos", "p.Neg"))) {
    stop("stats can only allow either cor.coef, t, p.Pos and p.Neg.")
  }
  if (all(cell.type != "all")) {
    cell.type <- janitor::make_clean_names(cell.type, case = "sentence")
    if (length(cell.type) < 2L) {
      stop("cell.type must be either all or length larger than 1.")
    }
    if (!all(cell.type %in% unique(c(fit_dat$celltype1, fit_dat$celltype2)))) {
      stop("cell.type must be all in celltype1 or celltype2 columns.")
    }
    fit_dat <- fit_dat[(fit_dat$celltype1 %in% cell.type) & 
                         (fit_dat$celltype2 %in% cell.type), ]
  }
  if (!("ROI" %in% colnames(fit_dat))) {
    fit_dat <- fit_dat
    title <- paste0("Statistics (", stats, ")")
    filled_data <- plot_type1_transformation(fit_dat, stats)
  }
  else {
    if (all(roi != "all")) {
      if (length(roi) == 1L) {
        if (!(roi %in% fit_dat[["ROI"]])) {
          stop("roi is not in the ROI column.")
        }
        fit_dat <- fit_dat[fit_dat$ROI == roi, ]
        title <- paste0("Statistics (", stats, ") of ROI #", 
                        roi)
        filled_data <- plot_type1_transformation(fit_dat, 
                                                 stats)
      }
      else {
        if (!all(roi %in% fit_dat[["ROI"]])) {
          stop("rois are not all in the ROI column.")
        }
        fit_dat <- fit_dat[fit_dat$ROI %in% roi, ]
        title <- paste0("Statistics (", stats, ") of ROI #", 
                        paste(roi, collapse = ","))
        filled_data <- plot_type2_transformation(fit_dat, 
                                                 stats)
      }
    }
    else {
      title <- paste0("Statistics (", stats, ")")
      filled_data <- plot_type2_transformation(fit_dat, 
                                               stats)
    }
  }
  paletteLength <- 100
  hmColor <- rev(colorRampPalette(c("#ED254EFF", "#EF6079FF", 
                                    "#F1F4FFFF", "#97B3D0FF", "#011936FF"))(paletteLength))
  myBreaks <- c(seq(min(filled_data), 0, length.out = ceiling(paletteLength/2) + 
                      1), seq(max(filled_data)/paletteLength, max(filled_data), 
                              length.out = floor(paletteLength/2)))
  tmp <- unique(paste0(model.result[,"celltype1"],model.result[,"celltype2"]))
  if (length(tmp)>200) {cellsize=5} else {
    if (length(tmp)>100) {cellsize=6}  else {cellsize=7} }
  if (nrow(filled_data) == 1L) {
    pheatmap::pheatmap(filled_data, angle_col = 45, border_color = "white", 
                       color = hmColor, breaks = myBreaks, main = title, 
                       cluster_rows = FALSE,fontsize_row = cellsize)
  }
  else {
    pheatmap::pheatmap(filled_data, angle_col = 45, border_color = "white", 
                       color = hmColor, breaks = myBreaks, main = title,fontsize_row = cellsize)
  }
}


preprocess<-function(SampleInfo, Cutoff,  Mtpattern="^MT-", resolution=0.5, Remove_Mt_rRNA=FALSE,celltype_predictmethod="cta",transpose=FALSE, analysis_individual=TRUE,Ensemblfile=NULL) {

  v3d <- Load10X_Spatial(
    data.dir = SampleInfo$countfile,
    slice = "slice1",
    bin.size = c(8, 16, "polygons")
  )
  v3d[['Spatial.Polygons.percent.mt']] <- PercentageFeatureSet(v3d,assay="Spatial.Polygons", pattern = Mtpattern)

  cat("\n\n# Sample", SampleInfo$SampleId,": Quality Check and Analysis\n\n")
  cat("\n\n## ", SampleInfo$SampleId,": Quality Check\n\n")
  cat("\n\n### ", "Fig.1 Violin plot of nGene,nUMI and mtRNA distribution\n\n")
  print(VlnPlot(v3d, features = c("nFeature_Spatial.Polygons","nCount_Spatial.Polygons","Spatial.Polygons.percent.mt"), ncol = 3, pt.size=0.01))
  cat("\n\n### ", "Fig.2 The scatterplot between mtRNA/nGene and nUMI\n\n")
  nCount_mt_plot <- FeatureScatter(v3d, feature1 = "nCount_Spatial.Polygons", feature2 = "Spatial.Polygons.percent.mt",pt.size=1) + ggtitle("The scatterplot between mtRNA and nUMI")
  nCount_nFeature_plot <- FeatureScatter(v3d, feature1 = "nCount_Spatial.Polygons", feature2 = "nFeature_Spatial.Polygons",pt.size=1) + ggtitle("The scatterplot between nGene and nUMI")
  print(CombinePlots(plots = list(nCount_mt_plot, nCount_nFeature_plot)))

  # Cutoff<-data.frame(nFeature_cutoff_min = 25 ,nFeature_cutoff_max = 1000,nCount_cutoff=50, mt_cutoff = 25, cluster_remove=c(""),stringsAsFactors = F)
  v3d <- subset(v3d, subset = nFeature_Spatial.Polygons > Cutoff$nFeature_cutoff_min &
                  nFeature_Spatial.Polygons < Cutoff$nFeature_cutoff_max &
                  nCount_Spatial.Polygons > Cutoff$nCount_cutoff &
                  Spatial.Polygons.percent.mt < Cutoff$mt_cutoff)

  cat("\n\n### ", "Fig.3 nGene distribution on the original tissue image\n\n")
  print(SpatialFeaturePlot(
    object = v3d,
    images = "slice1.polygons",
    features = "nFeature_Spatial.Polygons",
    plot_segmentations = TRUE,
    pt.size.factor = 0.2
  ))

  v3d <- SCTransform(v3d,assay = "Spatial.Polygons",new.assay.name = "Polygon")
  DefaultAssay(v3d) <- "Polygon"
  v3d <- RunPCA(v3d)
  v3d <- FindNeighbors(v3d, dims = 1:30)
  v3d <- FindClusters(v3d, resolution = resolution)
  v3d <- RunUMAP(v3d, reduction = "pca", dims = 1:30, resolution = resolution)
  v3d <- CellCycleScoring(v3d,s.genes,g2m.genes)
  # saveRDS(v3d,file=paste0(workdir,"/WHY_01.rds"))
  cat("\n\n### ", "Fig.4 nGene,nUMI and mtRNA distribution in each cluster\n\n")
  print(FeaturePlot(v3d, features=c("Spatial.Polygons.percent.mt","nFeature_Polygon","nCount_Polygon","PC_1"),label=T))
  cat("\n\n### ", "Fig.5 Boxplot of nUMI,nGene, mtRNA and nCells distribution and PCA, UMAP results\n\n")
  plot1<- ggplot(data=data.frame(cluster=v3d@active.ident,nUMI=v3d$nCount_Spatial.Polygons),aes(x=cluster,y=nUMI))+geom_boxplot()+scale_y_log10()
  plot2<- ggplot(data=data.frame(cluster=v3d@active.ident,ngene=v3d$nFeature_Spatial.Polygons),aes(x=cluster,y=ngene))+geom_boxplot()+scale_y_log10()
  plot3<-ggplot(data=data.frame(cluster=v3d@active.ident,MtRNA=v3d$Spatial.Polygons.percent.mt),aes(x=cluster,y=MtRNA))+geom_boxplot()
  plot4<- ggplot(data=data.frame(cluster=v3d@active.ident),aes(x=cluster))+geom_bar(stat="count")+theme(axis.text.x = element_text(angle = 45))+ylab("nCells")+ggtitle(paste0("total cells:",dim(v3d)[2]))
  plot5<- DimPlot(v3d, reduction = "pca",label=T,label.size=4)
  plot6<-DimPlot(v3d, reduction = "umap",label=T,label.size=4)
  print(CombinePlots(plots = list(plot1, plot2, plot3,plot4,plot5,plot6),rel_heights=c(1,1)))
  cat("\n\n## ", SampleInfo$SampleId,": Cell Cycle Scoring\n\n")
  cat("\n\n### ", "Fig.6 Cell Cycle plot\n\n")
  print(DimPlot(v3d,group.by="Phase"))
  cat("\n\n## ", SampleInfo$SampleId,": Annotations of cell clusters\n\n")
  medianexp<-t(apply(GetAssayData(v3d,slot="data"),1,function(x){tapply(x,v3d@active.ident,mean)}))
  if (nrow(medianexp)==1) medianexp<-t(medianexp)
  predict_celltype<-ORA_celltype(medianexp,cellType,method=celltype_predictmethod)
  new.celltype.ids<-paste(rownames(predict_celltype$predict_result),sep="_")
  new.cluster.ids<-paste(levels(v3d@active.ident),rownames(predict_celltype$predict_result),sep="_")
  names(new.celltype.ids) <- levels(v3d)
  names(new.cluster.ids) <- levels(v3d)
  v3d$celltype <- as.character(new.celltype.ids[as.character(v3d$seurat_clusters)])
  v3d <- RenameIdents(v3d, new.cluster.ids)
  v3d$seurat_clusters_ann <- Idents(v3d)
  cat("\n\n### ", "Fig.7 Cluster annotation\n\n")
  print(DimPlot(v3d, reduction = "umap",label=T,label.size=6))
  cat("\n\n### ", "Fig.8 cell type annotation\n\n")
  print(DimPlot(v3d, group.by="celltype",reduction = "umap",label=T,label.size=6))
  cat("\n\n### Fig. 9 cluster plot on original tissue \n\n")
  print(SpatialDimPlot(v3d, images = "slice1.polygons", pt.size.factor = 0.5, alpha = c(1, 1),image.alpha = 0,
                plot_segmentations = TRUE,
                facet.highlight = TRUE))
  allidents <- unique(Idents(v3d))
  allidents <- allidents[order(as.numeric(sub("_.*", "", allidents)))]
  for (j in 1:length(allidents)) {
    print(SpatialDimPlot(v3d, images = "slice1.polygons", pt.size.factor = 0.5, alpha = c(1, 1),image.alpha = 0,
                  cells.highlight = CellsByIdentities(object = v3d, idents = allidents[j]),
                  plot_segmentations = TRUE,
                  facet.highlight = TRUE))
  }
  cat("\n\n## ", SampleInfo$SampleId,": Markers of cell clusters\n\n")
  cat("\n\n### Fig. 10 Expression of user defined marker genes \n\n")
  print(get_bubble_plot(obj=v3d,
                        NULL,
                        NULL,
                        bubblemap_file=bubblemap_file,
                        assay="Polygon",
                        group.by = "seurat_clusters_ann",
                        species=species))
  v3d.markers <- FindAllMarkers(v3d, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  v3d@misc$markers<-v3d.markers
  top10 <- v3d.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  if (nrow(top10)>200) {genesize=5} else {
    if (nrow(top10)>100) {genesize=6}  else {genesize=7} }
  heatmap_of_markers1 <- DoHeatmap(v3d, features = top10$gene,size = 1) + 
    scale_fill_viridis_b() + theme(axis.text.y = element_text(size = genesize)) + NoLegend()
  ggsave(paste0(workdir,"/data/heatmap_of_markers_withinscaledata.png"),heatmap_of_markers1,width=2000, height=1600, dpi=300, units="px", bg="white")
  heatmap_of_markers2 <- DoHeatmap(v3d, features = top10$gene,slot="data",size = 1) + 
    scale_fill_viridis_b() + theme(axis.text.y = element_text(size = genesize)) + NoLegend()
  ggsave(paste0(workdir,"/data/heatmap_of_markers_withindata.png"),heatmap_of_markers2,width=2000, height=1600, dpi=300, units="px", bg="white")

  return(v3d)
}


#### density and colocalization


# scider_process <- function (SampleInfo,object) {
#   # object<-readRDS(file=paste0("data/",SampleInfo$SampleId,"_seurat.rds"))
#   spe <- seurat_to_spe(seu = object, sample_id = "SeuratProject", img_id = "slice1.polygons")
#   colData(spe)$cell_type <- colData(spe)$seurat_clusters_ann
#   cat("\n\n## ", SampleInfo$SampleId,": Density and colocalization analyses\n\n")
#   cat("\n\n### Fig. 12 Density of each cluster\n\n")
#   coi <- unique(colData(spe)$cell_type)
#   spe <- gridDensity(spe, coi = coi)
#   coi <- coi[order(as.numeric(sub("_.*", "", coi)))]
#   for (k in 1:length(coi)) {
#     print (plotDensity(spe,coi=coi[k])+ ggtitle(coi[k]))
#   }
#   cat("\n\n### Fig. 13 Regions-of-interst based on density\n\n")
#   spe <- findROI(spe, coi = coi)
#   print(plotROI(spe))
#   cat("\n\n### Fig. 14 Correlations between cell types in each Region-of-interst\n\n")
#   model_result <- corDensity(spe)
#   tmp <- plotCorHeatmap(model_result)
#   if (length(test$tree_row$labels)>200) {cellsize=4} else {
#     if (length(test$tree_row$labels)>100) {cellsize=5}  else {cellsize=6} }
#   print(plotCorHeatmap(model_result) + theme(axis.text.y = element_text(size = cellsize)))
#   cat("\n\n### Fig. 15 Cell annotation based on density\n\n")
#   for (j in 1:length(coi)) {
#     spe <- getContour(spe, coi = coi[j])
#     spe <- allocateCells(spe)
#     p1 <- plotContour(spe,coi = coi[j])
#     p2 <- plotCellCompo(spe, coi = coi[j])
#     print(CombinePlots(plots = list(p1, p2),nrow=2, rel_heights = c(3,1)))
#   }
#   return(spe)
# }


scider_process <- function (SampleInfo) {
  # object<-readRDS(file=paste0("data/",SampleInfo$SampleId,"_seurat.rds"))
  # spe <- seurat_to_spe(seu = object, sample_id = "SeuratProject", img_id = "slice1.polygons")
  spe <- readRDS(paste0("data/",SampleInfo$SampleId,"_spe.rds"))
  # colData(spe)$cell_type <- colData(spe)$seurat_clusters_ann
  # cat("\n\n## ", SampleInfo$SampleId,": Density and colocalization analyses\n\n")
  # cat("\n\n### Fig. 12 Density of each cluster\n\n")
  # coi <- unique(colData(spe)$cell_type)
  # spe <- gridDensity(spe, coi = coi)
  # coi <- coi[order(as.numeric(sub("_.*", "", coi)))]
  # for (k in 1:length(coi)) {
  #   print (plotDensity(spe,coi=coi[k])+ ggtitle(coi[k]))
  # }
  # cat("\n\n### Fig. 13 Regions-of-interst based on density\n\n")
  # spe <- findROI(spe, coi = coi)
  # print(plotROI(spe))
  cat("\n\n### Fig. 14 Correlations between cell types in each Region-of-interst\n\n")
  model_result <- corDensity(spe)
  print(plotCorHeatmap_modified(model.result=model_result))
  model_result2 <- corDensity(spe, by.roi = FALSE)
  print(plotCorHeatmap(model_result2))
  # cat("\n\n### Fig. 15 Cell annotation based on density\n\n")
  # for (j in 1:length(coi)) {
  #   spe <- getContour(spe, coi = coi[j])
  #   spe <- allocateCells(spe)
  #   p1 <- plotContour(spe,coi = coi[j])
  #   p2 <- plotCellCompo(spe, coi = coi[j])
  #   print(CombinePlots(plots = list(p1, p2),nrow=2, rel_heights = c(3,1)))
  # }
  # return(spe)
}






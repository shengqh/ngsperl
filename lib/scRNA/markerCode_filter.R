
library(dplyr)

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
Doheatmap_cellmarker_cellType<-function(SCLC,top10,cellType,predict_celltype) {
  genemarker<-rev(top10$gene)
  orderind<-order(SCLC@active.ident)
  data_m<-GetAssayData(SCLC, assay="RNA")
  geneind<-match(genemarker,rownames(data_m))
  data_h<-data_m[geneind,orderind]
  data_h_scale<-scale(t(as.matrix(data_h)))
  data_h_scale<-MinMax(data_h_scale,min=-2.5,max=2.5)
  
  matchid<-NULL
  ###
  for (i in 1:length(levels(SCLC@active.ident))){
    
    matchid<-c(matchid,ifelse(is.na(match(top10$gene[top10$cluster==i-1],cellType[[rownames(predict_celltype$predict_result)[i]]])),0,1))
    
  }
  matchid<-rev(matchid)
  predict_result<-data.frame(predict_celltype$predict_result)
  
  rowrange<-seq(from = 0, to = 1, length =ncol(data_h_scale))
  colrange<-seq(from = 0, to = 1, length =nrow(data_h_scale))
  clustnum<-table(SCLC@active.ident)
  colind<-c(0,cumsum(clustnum))+c(clustnum/2,0)
  colind<-floor(colind)[-length(colind)]
  
  col_text<-ifelse(predict_result$p_val<1e-5,"blue","gray")
  names(colind)<-levels(SCLC@active.ident)
  
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

get_mean_expression<-function(SCLC){
  result<-t(apply(GetAssayData(SCLC,slot="data"),1,function(x){tapply(x,SCLC@active.ident,mean)}))
  return(result)
}

get_seurat_average_expression<-function(SCLC){
  dd=GetAssayData(SCLC,slot="data")
  dobj=CreateSeuratObject(counts=dd)
  dobj$seurat_clusters=SCLC$seurat_clusters
  result<-AverageExpression(dobj, slot="counts", group.by="seurat_clusters" )[[1]]
  rm(dd)
  rm(dobj)
  return(result)
}

#https://github.com/yihui/knitr/issues/1494
# Based on code from http://michaeljw.com/blog/post/subchunkify/
#' Generate a sub-chunk to be interpreted by knitr.  The enclosing chunk
#' must have "results='asis'"
#'
#' @param g The output to chunkify (only tested with figures to date)
#' @param ... Additional named arguments to the chunk
#' @return NULL
#' @details The chunk is automatically output to the console.  There is
#'   no need to print/cat its result.
#' @export
subchunkify <- local({
  chunk_count <- 0
  function(g, ...) {
    chunk_count <<- chunk_count + 1
    g_deparsed <-
      paste0(deparse(
        function() {print(g)}
      ),
      collapse = '')
    args <- list(...)
    args <-
      lapply(names(args),
             FUN=function(nm, arglist) {
               current <- arglist[[nm]]
               if (length(current) > 1) {
                 stop("Only scalars are supported by subchunkify")
               } else if (is.character(current) | is.factor(current)) {
                 current <- as.character(current)
                 ret <- paste0('"', gsub('"', '\"', current, fixed=TRUE), '"')
               } else if (is.numeric(current) | is.logical(current)) {
                 ret <- as.character(current)
               } else {
                 stop("Unhandled class in subchunkify argument handling")
               }
               paste0(nm, "=", ret)
             },
             arglist=args)
    args <- paste0(unlist(args), collapse=", ")
    chunk_header <-
      paste(
        paste0("{r sub_chunk_", chunk_count),
        if (nchar(args) > 0) {
          paste(",", args)
        } else {
          NULL
        },
        ", echo=FALSE}")
    
    sub_chunk <- paste0(
      "\n```",chunk_header, "\n",
      "(", 
      g_deparsed
      , ")()\n",
      "```\n")
    cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
  }
})

myScaleData<-function(object, features, assay, ...){
  scaled.genes<-rownames(object[[assay]]@scale.data)
  if(!all(features %in% scaled.genes)){
    new.genes<-unique(features, scaled.genes)
    object=ScaleData(object, features=new.genes, assay=assay, ... )
  }
  return(object)
}

preprocess<-function(SampleInfo, Cutoff,  Mtpattern="^MT-", resolution=0.5, Remove_Mt_rRNA=FALSE,celltype_predictmethod="cta",transpose=FALSE, Ensemblfile=NULL, hto_map=list(),tag_tb=NULL,bubble_file="") {
  countfile<-SampleInfo$countfile
  sampleid=SampleInfo$SampleId
  
  has_bubble_file<-!is.null(bubble_file)
  if(has_bubble_file){
    has_bubble_file=file.exists(bubble_file)
  }
  
  library(patchwork)
  
  if(url.exists(countfile) | (file.exists(countfile) & !dir.exists(countfile))){
    if (length(grep (".rds$",countfile))>0) { 
      counts<-readRDS(countfile)$cm
    } else if (length(grep (".h5$",countfile))>0) {
      counts<-Read10X_h5(countfile)
      if (is.list(counts)){
        counts<-counts$`Gene Expression` 
      }
    } else  {
      counts = data.frame(fread(countfile),row.names = 1,check.names=F) 
    }
  } else {
    counts=Read10X(countfile)
  }
  if (transpose){
    counts<-t(counts)
  }
  
  if(!is.null(Ensemblfile)){
    counts<-convertEnsembltoGeneSymbol(counts,Ensemblfile)
  }
  
  hto_file=hto_map[[sampleid]]
  has_hto=!is.null(hto_file)
  if(has_hto){
    hto<-read.csv(hto_file, stringsAsFactors = F, row.names=1)   
    hto<-subset(hto, !(hto$HTO %in% c("Negative", "Doublet")))
    hto_tag=subset(tag_tb, tag_tb$File==sampleid)
    tag_map=split(hto_tag$Sample, hto_tag$Tagname)
    hto$sample<-unlist(tag_map[hto$HTO])
    counts<-counts[, colnames(counts) %in% rownames(hto)]
  }

  cat("\n\n#", sampleid, "\n\n")
  cat("\n\n## Quality Check\n\n")
  
  #counts<-counts[,sample(ncol(counts), 2000)]

  SCLC <- CreateSeuratObject(counts = counts, min.cells = 5, min.features = 10, project=sampleid)
  if(has_hto){
    samples<-hto[colnames(SCLC), "sample"]
  }else{
    samples=sampleid
  }
  SCLC$sample=samples
  rm(counts)

  SCLC[["percent.mt"]] <- PercentageFeatureSet(SCLC, pattern = Mtpattern)
  
  filters<-c(list(sample=sampleid), as.list(Cutoff))
  filters$raw_num_cell=ncol(SCLC)
  
  cat("\n\n### ", "Fig.1 Violin plot of nGene,nUMI and mtRNA distribution\n\n")
  if(has_hto){
    g1<-VlnPlot(SCLC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "sample")
  }else{
    g1<-VlnPlot(SCLC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  }

  subchunkify(g1, fig.height=7, fig.width=15)
  #print(g1)

  plot1 <- FeatureScatter(SCLC, feature1 = "nCount_RNA", feature2 = "percent.mt")+NoLegend()
  plot2 <- FeatureScatter(SCLC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+NoLegend()
  cat("\n\n### ", "Fig.2 The scatterplot between mtRNA/nGene and nUMI \n\n")
  p<-plot1+plot2
  subchunkify(p, fig.height=7, fig.width=14)
  #print(plot1+plot2)
  
  mt<-data.frame(mt=SCLC$percent.mt, Sample=SCLC$sample, nFeature=log10(SCLC$nFeature_RNA), nCount=log10(SCLC$nCount_RNA))
  plot3<-ggplot(mt, aes(x=nCount,y=mt) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    geom_vline(xintercept = log10(Cutoff$nCount_cutoff), color="red")  + 
    geom_hline(yintercept = Cutoff$mt_cutoff, color="red") +
    ylab("Percentage of mitochondrial") + xlab("log10(number of read)") +
    theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  if(has_hto){
    plot3<-plot3+facet_grid(~Sample)
  }
  
  plot4<-ggplot(mt, aes(x=nFeature,y=mt) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    geom_vline(xintercept = log10(Cutoff$nFeature_cutoff_min), color="red")  + 
    geom_vline(xintercept = log10(Cutoff$nFeature_cutoff_max), color="red")  + 
    geom_hline(yintercept = Cutoff$mt_cutoff, color="red") +
    ylab("Percentage of mitochondrial") + xlab("log10(number of gene)") +
    theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  if(has_hto){
    plot4<-plot4+facet_grid(~Sample)
  }
  
  cat("\n\n### ", "Fig.3  nUMI distribution\n\n")
  if(has_hto){
    subchunkify(plot3, fig.height=7, fig.width=15)
    subchunkify(plot4, fig.height=7, fig.width=15)
    #print(plot3)
    #print(plot4)
  }else{
    p<-plot3+plot4
    subchunkify(p, fig.height=7, fig.width=14)
  }


  # plot(sort(SCLC@meta.data$nCount_RNA,decreasing = T),ylab="UMI counts",lty=2,pch=16,log="xy")
  # abline(h=1000)
  # abline(h=500)
  
  info=list()
  if (min(SCLC[["percent.mt"]]) < Cutoff$mt_cutoff) {  
    SCLC <- subset(SCLC, subset = nFeature_RNA > Cutoff$nFeature_cutoff_min & nFeature_RNA<Cutoff$nFeature_cutoff_max & nCount_RNA>Cutoff$nCount_cutoff & percent.mt < Cutoff$mt_cutoff)
    filters$valid_num_cell=ncol(SCLC)

    SCLC <- NormalizeData(SCLC)
    SCLC <- FindVariableFeatures(SCLC, selection.method = "vst", nfeatures = 2000)
    
    rRNA.genes <- grep(pattern = rRNApattern,  rownames(SCLC), value = TRUE)
    Mt.genes<- grep (pattern= Mtpattern,rownames(SCLC), value=TRUE ) 
    
    if (Remove_Mt_rRNA) {
      var.genes <- dplyr::setdiff(VariableFeatures(SCLC), c(rRNA.genes,Mt.genes))
    } else {
      var.genes <- VariableFeatures(SCLC)
    }

    if (dim(SCLC)[2]>50){
      SCLC <- ScaleData(SCLC, vars.to.regress = c("percent.mt"))
      
      #find variable genes and store them in the var.genes
      SCLC <- RunPCA(SCLC, features = var.genes)
      SCLC <- FindNeighbors(SCLC, dims = 1:30)
      SCLC <- FindClusters(SCLC, resolution = resolution)
      ##
      SCLC[["celltype_each"]] <- paste(sampleid,SCLC@active.ident,sep="_")
      
      SCLC <- RunUMAP(SCLC, dims = 1:30)
      
      cat("\n\n## Cell clusters\n\n")
      cat("\n\n### ", "Fig.4 nGene,nUMI and mtRNA distribution in each cluster and PCA, UMAP results\n\n")
      if(has_hto){
        g1<-FeaturePlot(SCLC, features=c("percent.mt","nFeature_RNA","nCount_RNA","PC_1"),label=T, split.by = "sample")
      }else{
        g1<-FeaturePlot(SCLC, features=c("percent.mt","nFeature_RNA","nCount_RNA","PC_1"),label=T)
      }
      subchunkify(g1, fig.height=15, fig.width=15)
      
      plot1<-ggplot(data=data.frame(cluster=SCLC@active.ident,nUMI=SCLC@meta.data$nCount_RNA),aes(x=cluster,y=nUMI))+geom_boxplot()+scale_y_log10()    + theme_bw()
      plot2<-ggplot(data=data.frame(cluster=SCLC@active.ident,ngene=SCLC@meta.data$nFeature_RNA),aes(x=cluster,y=ngene))+geom_boxplot()+scale_y_log10() + theme_bw()
      plot3<-ggplot(data=data.frame(cluster=SCLC@active.ident,MtRNA=SCLC@meta.data$percent.mt),aes(x=cluster,y=MtRNA))+geom_boxplot()   + theme_bw()
      if(has_hto){
        plot4<-ggplot(data=data.frame(cluster=SCLC@active.ident, sample=SCLC$sample),aes(x=cluster, fill=sample))+geom_bar(stat="count")+theme(axis.text.x = element_text(angle = 45))+ylab("nCells")+ggtitle(paste0("total cells:",dim(SCLC)[2])) + theme_bw() + theme(legend.position="top")
      }else{
        plot4<-ggplot(data=data.frame(cluster=SCLC@active.ident),aes(x=cluster))+geom_bar(stat="count")+theme(axis.text.x = element_text(angle = 45))+ylab("nCells")+ggtitle(paste0("total cells:",dim(SCLC)[2])) + theme_bw()
      }
      plot5<-DimPlot(SCLC, reduction = "pca",label=T,label.size=4) 
      plot6<-DimPlot(SCLC, reduction = "umap",label=T,label.size=4)
      
      cat("\n\n### ", "Fig.5 Boxplot of nUMI,nGene, mtRNA and nCells distribution and PCA, UMAP results\n\n") 
      p<-plot1+plot2+plot3+plot4+plot5+plot6+plot_layout(ncol = 3)
      subchunkify(p, fig.height=10, fig.width=15)
      
      if (Cutoff$cluster_remove!=""){
        removeIdent<-(unlist((strsplit(Cutoff$cluster_remove,",")))) 
        SCLC<-subset(SCLC,idents=removeIdent, invert=T)
      }
      
      cat("\n\n## Markers and cell type annotation\n\n")

      ##predict method
      if (!is.null(celltype_predictmethod) & length(levels(SCLC@active.ident))>1 ) {
        meanexp<-get_seurat_average_expression(SCLC)
        if (nrow(meanexp)==1) meanexp<-t(meanexp)
        predict_celltype<-ORA_celltype(meanexp,cellType,method=celltype_predictmethod)
        
        #use the max_cta_score for each cluster
        new.cluster.ids<-paste(levels(SCLC@active.ident),rownames(predict_celltype$predict_result),sep="_")
        
        names(new.cluster.ids) <- levels(SCLC)
        SCLC <- RenameIdents(SCLC, new.cluster.ids)
        cat("\n\n### ", "Fig.6 UMAP result\n\n")
        g<-DimPlot(SCLC, reduction = "umap",label=T,label.size=6)
        subchunkify(g, fig.height=15, fig.width=15)

        if(has_hto){
          g<-DimPlot(SCLC, reduction = "umap",label=T,label.size=3, split.by = "sample")
          subchunkify(g, fig.height=15, fig.width=15)
        }
        
        SCLC.markers <- FindAllMarkers(SCLC, assay="RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
        SCLC@misc$markers<-SCLC.markers
        
        if('avg_log2FC' %in% colnames(SCLC.markers)){
          top10 <- SCLC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
        }else{
          top10 <- SCLC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
        }
        #cat("\n\n### ", "Fig.7 Marker genes scaled expression in each cluster\n\n")
        top10gene<-unique(top10$gene)
        if (length(top10gene)>200) {
          genesize=5
        } else {
          if (length(top10gene)>100) {
            genesize=6
          }  else {
            genesize=7
          } 
        }
        
        SCLC<-myScaleData(SCLC, top10gene, "RNA")
        
        #print(DoHeatmap(SCLC, features = top10$gene)+ theme(axis.text.y = element_text(size = genesize)) )
        cat("\n\n### Fig.7 Marker genes expression in each cluster\n\n")
        g<-DoHeatmap(SCLC, assay="RNA", features = top10gene)+ theme(axis.text.y = element_text(size = genesize))
        subchunkify(g, fig.height=15, fig.width=15)
        
        ##
        # SCLC@misc$Qcluster <- EvalCluster(SCLC)
        # ViewQcluster(SCLC)
        # cat("\n\n")
        # print(kable(SCLC@misc$Qcluster))
        # cat("\n\n")
        
        # Doheatmap_cellmarker_cellType(SCLC,top10,cellType,predict_celltype)
        # def.par <- par(no.readonly = TRUE) 
        #   layout(matrix(c(1,2), 1, 2, byrow = TRUE))
        #   Plot_predictcelltype(predict_celltype,method="cta")
        #   Plot_predictcelltype(predict_celltype,method="ora")
        # par(def.par)
        
        if(has_bubble_file){
          genes<-read_bubble_genes(bubble_file, rownames(SCLC))
          ugenes<-unique(genes$gene)
          
          cat("\n\n### Fig.8 Cell type marker genes expression in each cluster\n\n")
          
          SCLC<-myScaleData(SCLC, ugenes, "RNA")

          g<-DoHeatmap(SCLC, assay="RNA",features=ugenes)
          subchunkify(g, fig.height=15, fig.width=15)
          
          gene_groups=split(genes$gene, genes$cell_type)
          
          genes=unique(unlist(gene_groups))
          g<-DotPlot(SCLC, assay="RNA", features=genes)
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
          dot.scale = 6
          cols = c("lightgrey", "blue")
          
          library(cowplot)
          plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "id")) + 
            geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
            scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
            labs(x = "Features", y = "Identity") +
            theme_cowplot() + 
            facet_grid(facets = ~feature.groups, scales = "free_x", space = "free_x", switch = "y") + 
            theme(panel.spacing = unit(x = 1,units = "lines"), strip.background = element_blank()) + 
            scale_color_gradient(low = cols[1], high = cols[2])
          
          g=plot + 
            xlab("") + ylab("") + theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                                                     axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
                                                     strip.background = element_blank(),
                                                     strip.text.x = element_text(angle=90, hjust=0, vjust=0.5))
          cat("\n\n### Fig.9 Cell type marker genes bubble plot\n\n")
          subchunkify(g, fig.height=10, fig.width=15)
          cat("<br><br>\n\n\n")
        }
      }
    }

    info=list("preprocess"=filters, "meta"=SCLC@meta.data)
  }
  
  result=list(info)
  names(result)<-sampleid
  return(result)
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
  
  data<-GetAssayData(SCLC,slot="counts",assay="RNA")
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

getDifferential<-function(SumExp) {
  library(DESeq2)
  cluster<-levels(SumExp$group$cellcluster_ann)
  for (i in cluster) {
    tmp<-SumExp$Exp[,SumExp$group$cellcluster_ann==cluster[i]]
    tmp_group<-SumExp$group[SumExp$group$cellcluster_ann==cluster[i],]
    dds<-DESeq2::DESeqDataSetFromMatrix(countData=tmp,colData=tmp_group,design=~condition)
    dds <- DESeq(dds)
    res<-as.data.frame(results(dds))
    
    write.table(data.frame(res=res,celltype=cluster[i]),col.names = F,file="/Users/qiliu/Downloads/GSE155249_RAW (1)/analysis/diff_celltype.txt",sep="\t",quote=F,append = TRUE)
  }
}

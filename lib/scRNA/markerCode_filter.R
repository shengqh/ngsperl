
library(dplyr)
library(patchwork)

source("scRNA_func.r")

################

ORA_celltype_qc<-function(medianexp,cellType,method="cta"){
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

preprocess<-function( SampleInfo, 
                      Cutoff,  
                      cellType,
                      Mtpattern="^MT-|^Mt-|^mt-", 
                      rRNApattern="^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]",
                      resolution=0.5, 
                      Remove_Mt_rRNA=FALSE,
                      celltype_predictmethod="cta",
                      transpose=FALSE, 
                      Ensemblfile=NULL, 
                      hto_map=list(),
                      tag_tb=NULL,
                      bubblemap_file="",
                      bubblemap_width=6000,
                      bubblemap_height=3000) {

  countfile<-SampleInfo$countfile
  sampleid=SampleInfo$SampleId
  
  has_bubblemap_file<-!is.null(bubblemap_file)
  if(has_bubblemap_file){
    has_bubblemap_file=file.exists(bubblemap_file)
  }
  
  lst = read_scrna_data(countfile)
  counts<-lst$counts
  rm(lst)
  
  if (transpose){
    counts<-t(counts)
  }
  
  if(!is.null(Ensemblfile)){
    counts<-convertEnsembltoGeneSymbol(counts,Ensemblfile)
  }
  
  hto_file=hto_map[[sampleid]]
  has_hto=!is.null(hto_file)
  if(has_hto){
    if(grepl(".rds$", hto_file)){
      hto<-readRDS(hto_file)
    }else{
      hto<-read.csv(hto_file, stringsAsFactors = F, row.names=1)   
    }

    if("final" %in% colnames(hto)){
      colname="final"
    }else if ("HTO_classification" %in% colnames(hto)) {
      colname="HTO_classification"
    }else{
      stop("no final or HTO_classification in hto data, don't know which column should be used.")
    }

    hto<-subset(hto, !(hto[[colname]] %in% c("Negative", "Doublet")))
    hto_tag=subset(tag_tb, tag_tb$File==sampleid)
    tag_map=split(hto_tag$Sample, hto_tag$Tagname)
    hto$sample<-unlist(tag_map[unlist(hto[[colname]])])
    counts<-counts[, colnames(counts) %in% rownames(hto)]
  }

  cat("\n\n#", sampleid, "\n\n")
  #cat("\n\n## Quality Check\n\n")
  
  #counts<-counts[,sample(ncol(counts), 2000)]

  obj <- CreateSeuratObject(counts = counts, min.cells = 5, min.features = 10, project=sampleid)
  obj$sample=sampleid
  if(has_hto){
    obj$orig.ident=hto[colnames(obj), "sample"]
  }
  rm(counts)
  Idents(obj)<-"orig.ident"

  obj<-PercentageFeatureSet(object=obj, pattern=Mtpattern, col.name="percent.mt")
  obj<-PercentageFeatureSet(object=obj, pattern=rRNApattern, col.name = "percent.ribo")

  samples<-unique(obj$orig.ident)

  cur_folder = getwd()
  if(!dir.exists("details")){
    dir.create("details")
  }
  setwd(file.path(cur_folder, "details"))
  
  info=list()

  cur_sample=samples[1]
  for(cur_sample in samples){
    subobj=subset(obj, orig.ident==cur_sample)

    filters<-c(list(sample=cur_sample), as.list(Cutoff))
    filters$raw_num_cell=ncol(obj)
  
    #cat("\n\n### Violin plot of nGene,nUMI and mtRNA distribution\n\n")
    g1<-VlnPlot(subobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) & xlab("")
    png(paste0(cur_sample, ".qc1.png"), width=3000, height=1500, res=300)
    print(g1)
    dev.off()
  
    plot1 <- FeatureScatter(subobj, feature1 = "nCount_RNA", feature2 = "percent.mt")+NoLegend()
    plot2 <- FeatureScatter(subobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+NoLegend()
    #cat("\n\n### ", "Fig.2 The scatterplot between mtRNA/nGene and nUMI \n\n")
    p<-plot1+plot2
    png(paste0(cur_sample, ".qc2.png"), width=3000, height=1500, res=300)
    print(p)
    dev.off()
  
    mt<-data.frame(mt=subobj$percent.mt, nFeature=log10(subobj$nFeature_RNA), nCount=log10(subobj$nCount_RNA))
    plot3<-ggplot(mt, aes(x=nCount,y=mt) ) +
      geom_bin2d(bins = 70) + 
      scale_fill_continuous(type = "viridis") + 
      geom_vline(xintercept = log10(Cutoff$nCount_cutoff), color="red")  + 
      geom_hline(yintercept = Cutoff$mt_cutoff, color="red") +
      ylab("Percentage of mitochondrial") + xlab("log10(number of read)") +
      theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  
    plot4<-ggplot(mt, aes(x=nFeature,y=mt) ) +
      geom_bin2d(bins = 70) + 
      scale_fill_continuous(type = "viridis") + 
      geom_vline(xintercept = log10(Cutoff$nFeature_cutoff_min), color="red")  + 
      geom_vline(xintercept = log10(Cutoff$nFeature_cutoff_max), color="red")  + 
      geom_hline(yintercept = Cutoff$mt_cutoff, color="red") +
      ylab("Percentage of mitochondrial") + xlab("log10(number of gene)") +
      theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  
    #cat("\n\n### ", "Fig.3  nUMI distribution\n\n")
    p<-plot3+plot4
    png(paste0(cur_sample, ".qc3.png"), width=3000, height=1200, res=300)
    print(p)
    dev.off()

    subobj <- subset(subobj, subset = nFeature_RNA > Cutoff$nFeature_cutoff_min & nFeature_RNA<Cutoff$nFeature_cutoff_max & nCount_RNA>Cutoff$nCount_cutoff & percent.mt < Cutoff$mt_cutoff)
    filters$valid_num_cell=ncol(subobj)

    if(ncol(subobj) == 0){
      writeLines(paste0(cur_sample, " has no cells after filter."), paste0(cur_sample, ".nocell.txt"))
      next
    }

    subobj <- NormalizeData(subobj)
    subobj <- FindVariableFeatures(subobj, selection.method = "vst", nfeatures = 2000)
    
    rRNA.genes <- grep(pattern = rRNApattern,  rownames(subobj), value = TRUE)
    Mt.genes<- grep (pattern= Mtpattern,rownames(subobj), value=TRUE ) 
    
    if (Remove_Mt_rRNA) {
      var.genes <- dplyr::setdiff(VariableFeatures(subobj), c(rRNA.genes,Mt.genes))
    } else {
      var.genes <- VariableFeatures(subobj)
    }

    subobj <- ScaleData(subobj, vars.to.regress = c("percent.mt"))
    
    #find variable genes and store them in the var.genes
    subobj <- RunPCA(subobj, features = var.genes)
    subobj <- FindNeighbors(subobj, dims = 1:20)
    subobj <- FindClusters(subobj, resolution = resolution)
    subobj <- RunUMAP(subobj, dims = 1:20)
    
    #cat("\n\n## Cell clusters\n\n")
    #cat("\n\n### ", "Fig.4 nGene,nUMI and mtRNA distribution in each cluster and PCA, UMAP results\n\n")
    g1<-FeaturePlot(subobj, features=c("percent.mt","nFeature_RNA","nCount_RNA","PC_1"),label=T)
    png(paste0(cur_sample, ".qc4.png"), width=3000, height=2600, res=300)
    print(g1)
    dev.off()

    plot1<-ggplot(data=data.frame(cluster=subobj@active.ident,nUMI=subobj@meta.data$nCount_RNA),aes(x=cluster,y=nUMI))+geom_boxplot()+scale_y_log10()    + theme_bw()
    plot2<-ggplot(data=data.frame(cluster=subobj@active.ident,ngene=subobj@meta.data$nFeature_RNA),aes(x=cluster,y=ngene))+geom_boxplot()+scale_y_log10() + theme_bw()
    plot3<-ggplot(data=data.frame(cluster=subobj@active.ident,MtRNA=subobj@meta.data$percent.mt),aes(x=cluster,y=MtRNA))+geom_boxplot()   + theme_bw()
    plot4<-ggplot(data=data.frame(cluster=subobj@active.ident),aes(x=cluster))+geom_bar(stat="count")+theme(axis.text.x = element_text(angle = 45))+ylab("nCells")+ggtitle(paste0("total cells:",dim(subobj)[2])) + theme_bw()

    p<-plot1+plot2+plot3+plot4+plot_layout(ncol=4)
    png(paste0(cur_sample, ".qc5.png"), width=4000, height=1000, res=300)
    print(p)
    dev.off()

    plot5<-DimPlot(subobj, reduction = "pca",label=T,label.size=4) 
    plot6<-DimPlot(subobj, reduction = "umap",label=T,label.size=4)
    
    #cat("\n\n### ", "Fig.5 Boxplot of nUMI,nGene, mtRNA and nCells distribution and PCA, UMAP results\n\n") 
    p<-plot5+plot6+plot_layout(ncol = 2)
    png(paste0(cur_sample, ".qc6.png"), width=2200, height=1000, res=300)
    print(p)
    dev.off()

    if (Cutoff$cluster_remove!=""){
      removeIdent<-(unlist((strsplit(Cutoff$cluster_remove,",")))) 
      subobj<-subset(subobj,idents=removeIdent, invert=T)
    }
    
    #cat("\n\n## Markers and cell type annotation\n\n")

    ##predict method
    if (!is.null(celltype_predictmethod) & length(levels(subobj@active.ident))>1 ) {
      meanexp<-get_seurat_average_expression(subobj, "seurat_clusters")
      if (nrow(meanexp)==1) meanexp<-t(meanexp)
      predict_celltype<-ORA_celltype_qc(meanexp,cellType,method=celltype_predictmethod)
      saveRDS(predict_celltype, paste0(cur_sample, ".ct.rds"))
      
      #use the max_cta_score for each cluster
      new.cluster.ids<-rownames(predict_celltype$predict_result)
      names(new.cluster.ids) <- levels(subobj)
      subobj$cell_type=new.cluster.ids[subobj$seurat_clusters]
      
      new.cluster.ids<-paste0(levels(subobj$seurat_clusters), ": ", rownames(predict_celltype$predict_result))
      names(new.cluster.ids) <- levels(subobj)
      subobj$seurat_cell_type=factor(new.cluster.ids[subobj$seurat_clusters], levels=new.cluster.ids)

      Idents(subobj)<-"seurat_cell_type"
      
      #cat("\n\n### ", "Fig.6 UMAP result\n\n")
      g<-DimPlot(subobj, reduction = "umap",label=T,label.size=4) + NoLegend()
      png(paste0(cur_sample, ".qc7.png"), width=1500, height=1500, res=300)
      print(g)
      dev.off()

      subobj.markers <- FindAllMarkers(subobj, assay="RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
      subobj@misc$markers<-subobj.markers
      
      if('avg_log2FC' %in% colnames(subobj.markers)){
        top10 <- subobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
      }else{
        top10 <- subobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
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
      
      subobj<-myScaleData(subobj, top10gene, "RNA")
      
      #print(DoHeatmap(SCLC, features = top10$gene)+ theme(axis.text.y = element_text(size = genesize)) )
      #cat("\n\n### Fig.7 Marker genes expression in each cluster\n\n")
      g<-DoHeatmap(subobj, assay="RNA", features = top10gene)+ theme(axis.text.y = element_text(size = genesize))
      png(paste0(cur_sample, ".heatmap.png"), width=4000, height=4000, res=300)
      print(g)
      dev.off()

      # subobj@misc$Qcluster <- EvalCluster(subobj)
      # ViewQcluster(SCLC)
      # cat("\n\n")
      # print(kable(subobj@misc$Qcluster))
      # cat("\n\n")
      
      # Doheatmap_cellmarker_cellType(SCLC,top10,cellType,predict_celltype)
      # def.par <- par(no.readonly = TRUE) 
      #   layout(matrix(c(1,2), 1, 2, byrow = TRUE))
      #   Plot_predictcelltype(predict_celltype,method="cta")
      #   Plot_predictcelltype(predict_celltype,method="ora")
      # par(def.par)
      
      if(has_bubblemap_file){
        genes<-read_bubble_genes(bubblemap_file, rownames(subobj))
        ugenes<-unique(genes$gene)
        
        #cat("\n\n### Fig.8 Cell type marker genes expression in each cluster\n\n")
        
        subobj<-myScaleData(subobj, ugenes, "RNA")

        g<-DoHeatmap(subobj, assay="RNA",features=ugenes)
        png(paste0(cur_sample, ".bubble_heatmap.png"), width=4000, height=4000, res=300)
        print(g)
        dev.off()

        g<-get_bubble_plot(subobj, NULL, NULL, bubblemap_file, assay="RNA", group.by="seurat_cell_type") + theme(text = element_text(size=20))
        #cat("\n\n### Fig.9 Cell type marker genes bubble plot\n\n")
        png(paste0(cur_sample, ".bubble.png"), width=bubblemap_width, height=bubblemap_height, res=300)
        print(g)
        dev.off()
      }
    }

    info[[cur_sample]]=list("preprocess"=filters, "meta"=subobj@meta.data)
  }
  
  setwd(cur_folder)

  return(info)
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

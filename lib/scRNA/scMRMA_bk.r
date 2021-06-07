library(Seurat)
library(scMRMA)
library(plyr)

ORA_celltype<-function(medianexp,cellType,weight){
  ORA_result<-matrix(NA, nrow=length(cellType),ncol=dim(medianexp)[2])
  CTA_result<-matrix(0,nrow=length(cellType),ncol=dim(medianexp)[2])
  exp_z<-scale(medianexp)
  genenames<-rownames(medianexp)
  for (j in 1: dim(medianexp)[2]){
    clusterexp<-medianexp[,j]
    clusterexp_z<-exp_z[,j]
    for (i in 1:length(cellType)){

      ct_exp<-length(intersect(genenames[clusterexp>0],cellType[[i]]))
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
      CTA_result[i,j]<-sum(exp_ss*weight_ss)/(length(exp_ss)^(1/3))
    }
  }
  rownames(ORA_result)<-rownames(CTA_result)<-names(cellType)
  minp_ora_ind<- apply(ORA_result,2,function(x){which.min(x)})
  minp_ora<-apply(ORA_result,2,min)
  names(minp_ora)<-rownames(ORA_result)[minp_ora_ind]

  max_cta_ind<- apply(CTA_result,2,function(x){which.max(x)})
  max_cta<-apply(CTA_result,2,max,na.rm=T)
  names(max_cta)<-rownames(CTA_result)[max_cta_ind]
  return(list(ora=ORA_result,cta=CTA_result,min_ora=minp_ora,max_cta=max_cta))
}

get_celltype_weight<-function(marker, i) {
  celltype <- tapply(marker[,1],marker[,i],list)
  celltype <- lapply(celltype, unique)
  freq<-sort((table(unlist(celltype)))/length(celltype))
  if(max(freq)-min(freq)==0){
    weight <- rep(1,length(freq))
    names(weight) <- names(freq)
  }else{
    weight<-1+sqrt((max(freq)-freq)/(max(freq)-min(freq)))
  }
  return(list(celltype=celltype, weight=weight))
}

#get annotation result for each cluster
ori_annotation_normdata <- function(marker,i,data.norm,p){
  ct<-get_celltype_weight(marker, i)
  result<-ORA_celltype(data.norm,ct$celltype,ct$weight)
  colnames(result$ora)<-colnames(data.norm)
  colnames(result$cta)<-colnames(data.norm)
  result$ct=ct
  return(result)
}

get_annotation_list <-function(predict_celltype) {
  annoP <- predict_celltype$ora
  max_cta_pvalue = unlist(lapply(c(1:length(predict_celltype$max_cta)), function(x){
    cluster = x
    cta = names(predict_celltype$max_cta)[x]
    annoP[cta, cluster]
  }))
  celltype = unlist(lapply(c(1:length(max_cta_pvalue)), function(x){
    ifelse(max_cta_pvalue[x] > p, "unassigned", names(predict_celltype$max_cta)[x])
  }))

  result = data.frame(  predict_celltype=names(predict_celltype$max_cta),
                        activity= predict_celltype$max_cta,
                        pvalue=max_cta_pvalue,
                        celltype=celltype)
  rownames(result)<-colnames(predict_celltype$cta)

  return(result)
}

get_annotation_normdata <- function(data.norm, species,db="panglaodb",p=0.05){
  if(is(db,"data.frame")){
    CellType <- db
  }else if (is(db,"matrix")){
    CellType <- data.frame(db)
  }else {
    CellType <- scMRMA:::get_celltype(species,db)
  }

  annoResult <- list()
  anno <- matrix("unknown",nrow = dim(data.norm)[2],ncol = (ncol(CellType)) * 3)
  rownames(anno)<- colnames(data.norm)
  colnames(anno) <- paste(rep(c(colnames(CellType)[2:ncol(CellType)], "extend"),each=3),
                          rep(c("celltype","activity","pValue"),3),sep = "_")

  annoDetails <- matrix(NA,nrow = dim(data.norm)[2],ncol = (ncol(CellType)-1)*2)
  rownames(annoDetails)<- colnames(data.norm)
  colnames(annoDetails)<- paste(rep(colnames(CellType)[2:ncol(CellType)],each=2),
                                rep(c("celltypeActivity","pValue"),2),sep = "_")
  annoDetails <- as.data.frame(annoDetails,stringsAsFactors = FALSE)

  anno_markers<-list()
  anno<-list()
  cat("Multi Resolution Annotation Started.","\n")
  i=2
  lastlayer_data = NULL
  for (i in 2:ncol(CellType)){
    layer=paste0("layer", (i-1))
    cat(layer,"annotation started.","\n")
    if(i==2){
      marker <- CellType
      predict_celltype <- ori_annotation_normdata(marker,i,data.norm,p)
      lastlayer_anno = get_annotation_list(predict_celltype)
      anno[[layer]] = lastlayer_anno
    }else{
      curlayer_anno=NULL
      m= "CD4"
      for (m in unique(lastlayer_anno$celltype)){
        if (m != "unassigned"){
          cat("sub annotation of ", m, "\n")
          marker <- subset(CellType,CellType[,i-1]==m)
          if(length(unique(marker[,i])) > 1 ){##stop if there is no subtype
            clusternames <- rownames(lastlayer_anno)[lastlayer_anno$celltype==m]
            sub <- data.norm[,clusternames,drop=F]
            predict_celltype <- ori_annotation_normdata(marker,i,sub,p)
            sub_anno = get_annotation_list(predict_celltype)
            curlayer_anno=rbind(curlayer_anno, sub_anno)
          }else{
            sub_anno = subset(lastlayer_anno, lastlayer_anno$celltype == m)
            curlayer_anno=rbind(curlayer_anno, sub_anno)
          }
        }else{
          sub_anno = subset(lastlayer_anno, lastlayer_anno$celltype == m)
          curlayer_anno=rbind(curlayer_anno, sub_anno)
        }
      }
      anno[[layer]] = curlayer_anno
      lastlayer_anno = curlayer_anno
    }
  }

  predict_celltype <- ori_annotation_normdata(CellType,ncol(CellType),data.norm,p)
  anno[["uniformR"]] = get_annotation_list(predict_celltype)

  return(anno)
}

finalList<-readRDS(parFile1)
obj<-finalList$obj

data.norm<-read.csv(parFile2, header=T, row.names=1, check.names=F)
colnames(data.norm)<-gsub("Cluster","",colnames(data.norm))

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
options<-split(options_table$V1, options_table$V2)

species=options$species
prefix<-options$prefix
db=options$db
p=as.numeric(options$p)

res=get_annotation_normdata(data.norm, species, db, p)

for(layer in names(res)){
  curlayer = res[[layer]]
  curcelltype=plyr::mapvalues(obj$seurat_clusters, rownames(curlayer), as.character(curlayer$celltype))

  layer_name = paste0("scMRMA_", layer)
  obj[[layer_name]] <- curcelltype

  png(paste0(prefix, layer_name, ".png"), width=4000, height=3000, res=300)
  print(DimPlot(obj,reduction = "umap",group.by = layer_name,label = TRUE,repel = TRUE)+ggplot2::ggtitle(layer_name))
  dev.off()
}

finalListFile=paste0(prefix, ".scMRMA.rds")
saveRDS(finalList, file=finalListFile)

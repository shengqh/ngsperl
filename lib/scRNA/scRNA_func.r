read_cell_markers_file<-function(panglao5_file, species, remove_subtype_of="", HLA_panglao5_file=""){
  #preparing cell activity database
  marker<-data.frame(fread(panglao5_file))
  if(remove_subtype_of != ""){
    remove_subtype_of<-unlist(strsplit(remove_subtype_of, ","))
    pangdb_ct <- read.table(HLA_panglao5_file,header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
    layer="Layer2"
    removed<-c()
    for(layer in c("Layer1","Layer2", "Layer3")){
      layer_ct<-unique(pangdb_ct[,layer])
      for(rs in remove_subtype_of){
        if(rs %in% removed){
          next
        }
        
        if(rs %in% layer_ct){
          subdb<-rownames(pangdb_ct)[pangdb_ct[,layer]==rs]
          if(rs %in% subdb){ 
            subdb<-subdb[subdb != rs]
            marker<-marker[!(marker$cell.type %in% subdb),]
          }else{
            marker$cell.type[marker$cell.type %in% subdb] = rs
          }
          removed<-c(removed, rs)
        }
      }
    }
  }
  
  hsind<-regexpr(species,marker$species)
  marker_species<-marker[hsind>0 & marker$ubiquitousness.index<0.05,]
  if (species=="Mm") {
    ##change the gene symbol only keep the first letter capitalize
    marker_species$official.gene.symbol<-toMouseGeneSymbol(marker_species$official.gene.symbol)
  }
  if (species=="Hs") {
    marker_species$official.gene.symbol<-toupper(marker_species$official.gene.symbol)
  }
  cellType<-tapply(marker_species$official.gene.symbol,marker_species$cell.type,list)
  freq<-sort((table(unlist(cellType)))/length(cellType))
  weight<-1+sqrt((max(freq)-freq)/(max(freq)-min(freq)))
  return(list(cellType=cellType, weight=weight))
}

toMouseGeneSymbol<-function(x){
  result=paste0(toupper(substr(x,1,1)),tolower(substr(x,2,nchar(x))))
  return(result)
}

read_cell_cluster_file<-function(fileName, sort_cluster_name="seurat_clusters"){
  result<-read.csv(fileName, stringsAsFactors = F, row.names = 1)
  
  display_sort_cluster_name = paste0("display_", sort_cluster_name)
  result[,display_sort_cluster_name] = paste0("Cluster ", result[,sort_cluster_name])
  
  cluster_names=colnames(result)[grepl("_clusters", colnames(result))]
  
  sort_clusters_num = length(unique(result[,sort_cluster_name]))
  for(cluster_name in cluster_names){
    cluster_num = length(unique(result[,cluster_name]))
    if(cluster_name == sort_cluster_name){
      next
    }

    if (cluster_num != sort_clusters_num) {
      next
    }
      
    cf<-unique(result[,c(sort_cluster_name, cluster_name)])
    if(nrow(cf) != sort_clusters_num){
      next
    }
    
    cf<-cf[order(as.numeric(cf[,sort_cluster_name]), decreasing = T),]
    cf_levels=cf[,cluster_name]
    result[,cluster_name] = factor(result[,cluster_name], levels=cf_levels)
  }
  return(result)
}

find_markers<-function(object, by_sctransform, ident.1, ident.2, min.pct = 0.5, logfc.threshold = 0.6){
  assay=ifelse(by_sctransform, "SCT", "RNA")
  markers=FindMarkers(object, assay=assay, ident.1=ident.1, ident.2=ident.2, only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
  markers=markers[markers$p_val_adj < 0.01,]
  return(markers)
}

find_all_markers<-function(object, by_sctransform, min.pct = 0.5, logfc.threshold = 0.6){
  assay=ifelse(by_sctransform, "SCT", "RNA")
  markers=FindAllMarkers(object, assay=assay, only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
  markers=markers[markers$p_val_adj < 0.01,]
  return(markers)
}

get_cluster_count<-function(counts, clusters){
  if(is.null(levels(clusters))){
    allClusters=unique(clusters)
  }else{
    allClusters=levels(clusters)
  }
  cluster<-allClusters[1]
  csums=lapply(allClusters, function(cluster){
    #cat(cluster, "\n")
    cells=names(clusters)[clusters==cluster]
    subcounts=counts[,cells]
    #cat(ncol(subcounts), "\n")
    Matrix::rowSums(subcounts)
  })
  
  result=do.call(cbind, csums)
  colnames(result)=allClusters
  
  gcount=Matrix::rowSums(result)
  result=result[gcount > 0,]
  return(result)
}

get_group_count=function(curobj, groupName="active.ident") {
  counts=GetAssayData(curobj,assay="RNA",slot="counts")
  curgroups=curobj[[groupName]]
  clusters=curgroups[,1]
  names(clusters)=rownames(curgroups)
  result=get_cluster_count(counts, clusters)
  return(result)
}

add_cluster<-function(object, new.cluster.name, new.cluster.ids){
  seurat_clusters<-object[["seurat_clusters"]]$seurat_clusters
  names(new.cluster.ids) <- levels(seurat_clusters)

  new.cluster.values<-plyr::mapvalues(x = seurat_clusters, from = levels(seurat_clusters), to = new.cluster.ids)
  names(new.cluster.values)<-names(seurat_clusters)
  
  object[[new.cluster.name]]<-new.cluster.values
  object
}

run_cluster<-function(object, Remove_Mt_rRNA, rRNApattern, Mtpattern, pca_dims, by_sctransform, min.pct = 0.5, logfc.threshold = 0.6){
  if (by_sctransform) {
    object <- RunPCA(object = object, verbose=FALSE)
  }else{
    if (Remove_Mt_rRNA) {
      rRNA.genes <- grep(pattern = rRNApattern,  rownames(object), value = TRUE)
      Mt.genes<- grep (pattern= Mtpattern,rownames(object), value=TRUE )
      var.genes <- dplyr::setdiff(VariableFeatures(object), c(rRNA.genes,Mt.genes))
    } else {
      var.genes <- VariableFeatures(object)
    }
    object <- RunPCA(object = object, features = var.genes, verbose=FALSE)
  }
  object <- RunUMAP(object = object, dims=pca_dims, verbose = FALSE)
  object <- FindNeighbors(object = object, dims=pca_dims, verbose=FALSE)
  object <- FindClusters(object=object, verbose=FALSE, random.seed=random.seed, resolution=resolution)

  if (by_sctransform) {
    markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)
  }else{
    markers <- FindAllMarkers(object, features=var.genes, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)
  }
  markers <- markers[markers$p_val_adj < 0.01,]
  return(list(object=object, markers=markers))
}

ORA_celltype<-function(medianexp,cellType,weight){
  ORA_result<-matrix(NA, nrow=length(cellType),ncol=dim(medianexp)[2])
  CTA_result<-matrix(0,nrow=length(cellType),ncol=dim(medianexp)[2])

  colnames(ORA_result)=colnames(medianexp)
  colnames(CTA_result)=colnames(medianexp)

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

get_cta_ora_mat<-function(predict_celltype) {
  cta_index<-apply(predict_celltype$cta,2,function(x){return(order(x,decreasing=T)[1:2])})
  cta_index<-unique(sort(cta_index))
  
  cta_mat<- predict_celltype$cta[cta_index,]
  colnames(cta_mat)<-paste0(names(predict_celltype$max_cta), " : ", colnames(predict_celltype$cta))
  
  ora_mat<- predict_celltype$ora[cta_index,]
  ora_mat<--log10(ora_mat)  
  colnames(ora_mat)<-colnames(cta_mat)
  
  return(list(cta_mat=cta_mat, ora_mat=ora_mat))
}

get_cta_combined<-function(obj, predicted){
  cta_table<-data.frame(Cluster=colnames(predicted$cta), 
                        CellType=names(predicted$max_cta),
                        CtaScore=round(predicted$max_cta * 10) / 10.0,
                        stringsAsFactors = F)
  cta_table$OraPvalue=apply(cta_table, 1, function(x){
    ct=x[2]
    cl=x[1]
    predicted$ora[ct, cl]
  })

  cluster_sample<-as.data.frame.matrix(table(obj$seurat_clusters, obj$orig.ident))
  cluster_sample<-cluster_sample[as.character(cta_table$Cluster),]

  nc<-apply(cluster_sample, 2, function(x){
    tc=sum(x)
    perc<-x/tc
    return(round(perc*1000) / 10.0)
  })
  colnames(nc)<-paste0(colnames(nc), "_perc")

  cta_combined<-cbind(cta_table, cluster_sample, nc)

  return(cta_combined)
}

get_selfdefinedCelltype <- function(file){
  file.ori <- scan(file,what = "",sep="\n")
  marker.ori <- matrix("undefined",nrow=length(file.ori)/2,ncol = 2)
  for (i in 1:length(file.ori)){if(i%%2==1){marker.ori[(i+1)/2,1]=substr(file.ori[i],2,nchar(file.ori[i]))}else{marker.ori[i/2,2]=file.ori[i]}}
  ref <- as.data.frame(marker.ori[,1])
  ref <- tidyr::separate(ref,col=colnames(ref),into = c("celltype","subtypeOf"),sep = ",")
  rownames(ref) <- ref$celltype
  celltype.ori <- sapply(tapply(marker.ori[,2],marker.ori[,1],list), function(x) strsplit(x,","))

  celltype.tem <- as.data.frame(matrix(unlist(celltype.ori)))
  celltype.tem[,2] <- rep(names(celltype.ori),lengths(celltype.ori))
  celltype.tem <- tidyr::separate(celltype.tem,col=V2,into = c("celltype","subtypeOf"),sep = ",")

  ct <- as.data.frame(matrix("Undefined",nrow = nrow(celltype.tem),ncol = length(unique(celltype.tem$subtypeOf))+2),stringsAsFactors=FALSE)
  ct[,1:2] <- celltype.tem[,1:2]
  for (i in 3:ncol(ct)){
    for (j in 1:nrow(ct)) {if(ct[j,i-1] %in% rownames(ref)){ct[j,i]<-ref[ct[j,i-1],2]}else{ct[j,i]<-ct[j,i-1]}}
  }

  layer <- ncol(ct)
  if(! identical(ct[,ncol(ct)],ct[,ncol(ct)-1])){layer <- ncol(ct)}else{
    for (i in 1:ncol(ct)) {if(identical(ct[,i],ct[,i+1])){layer <- i;break}}
  }

  celltype.trim <- ct[,1:layer]
  for (i in 1:nrow(celltype.trim)){
    tag <- which(celltype.trim[i,]==celltype.trim[i,ncol(celltype.trim)])
    if(length(tag)>1){
      if(ncol(celltype.trim)+1-length(tag)==3){suf <- ct[i,3]}
      else{suf <- ct[i,3:(ncol(celltype.trim)+1-length(tag))]}
      pre <- rep(celltype.trim[i,2],length(tag)-1)
      celltype.trim[i,] <- c(celltype.trim[i,][1:2],pre,suf)
    }
  }

  for (i in 1:ncol(celltype.trim)){if(i == 1){colnames(celltype.trim)[i]="gene"}else{colnames(celltype.trim)[i]=paste0("layer",(ncol(celltype.trim)+1-i),"")}}

  ct.final <- celltype.trim
  for (i in 1:ncol(ct.final)){
    if(i >1){ct.final[,i]=celltype.trim[,ncol(celltype.trim)+2-i];colnames(ct.final)[i] <- colnames(celltype.trim)[ncol(celltype.trim)+2-i]}
  }

  cellType<-tapply(ct.final$gene,ct.final$layer3,list)
  freq<-sort((table(unlist(cellType)))/length(cellType))
  weight<-1+sqrt((max(freq)-freq)/(max(freq)-min(freq)))
  return(list(cellType=cellType, weight=weight))
}


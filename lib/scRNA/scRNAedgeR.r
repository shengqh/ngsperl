
library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)

sumcount<-function(ct_count, names, sample_df){
  result<-lapply(names, function(x){
    res<-ct_count[,sample_df$Cell[sample_df$Sample ==x],drop=F]
    apply(res, 1, sum)
  })
  rescount<-do.call(cbind, result)
  colnames(rescount)<-names
  return(rescount)
}

finalList<-readRDS(parFile1)
obj<-finalList$obj

clusterDf<-read.csv(parFile2, stringsAsFactors = F, row.names=1)
if(!(cluster_name %in% colnames(obj@meta.data))){
  if(all(names(obj$orig.ident) %in% rownames(clusterDf))){
    obj[[cluster_name]]<-clusterDf[names(obj$orig.ident), cluster_name]
  }else{
    obj[[cluster_name]]<-clusterDf[obj$seurat_clusters, cluster_name]
  }
}

meta<-obj@meta.data

comparisons<-read.table(parSampleFile2, stringsAsFactors = F)
if(ncol(comparisons) == 3){
  colnames(comparisons)<-c("Value", "Key", "Comparison")
}else{
  colnames(comparisons)<-c("Value", "Comparison")
  comparisons$Key = "groups"
}

comparisonNames<-unique(comparisons$Comparison)

comp <-comparisonNames[1]

designMatrix<-NULL
for (comp in comparisonNames){
  comp_groups<-comparisons[comparisons$Comparison==comp,]
  comp_options = split(comp_groups$Value, comp_groups$Key)
  
  if("groups" %in% names(comp_options)){
    sampleGroups<-read.table(parSampleFile1, stringsAsFactors = F)
    colnames(sampleGroups)<-c("Sample","Group")
    groups<-comp_options$groups
    controlGroup<-groups[1]
    sampleGroup<-groups[2]
    
    control_names<-sampleGroups$Sample[sampleGroups$Group==controlGroup]
    sample_names<-sampleGroups$Sample[sampleGroups$Group==sampleGroup]
  }else{
    control_names<-as.numeric(comp_options$control_clusters)
    controlGroup<-ifelse("control_name" %in% names(comp_options), comp_options$control_name, paste("Cluster", paste(control_names, collapse = "_"), sep="_"))
    sample_names<-as.numeric(comp_options$sample_clusters)
    sampleGroup<-ifelse("sample_name" %in% names(comp_options), comp_options$sample_name, paste("Cluster", paste(sample_names, collapse = "_"), sep="_"))
  }

  if(bBetweenCluster){
    prefix<-comp
    
    control_cells<-rownames(meta)[meta[,cluster_name] %in% control_names]  
    if("control_file_regex" %in% names(comp_options)){
      all_files = unique(meta$orig.ident)
      control_files=all_files[grepl(comp_options$control_file_regex, all_files)]
      control_cells = control_cells[meta[control_cells, "orig.ident"] %in% control_files]
    }
    if("control_files" %in% names(comp_options)){
      control_files=comp_options$control_files
      control_cells = control_cells[meta[control_cells, "orig.ident"] %in% control_files]
    }

    sample_cells<-rownames(meta)[meta[,cluster_name] %in% sample_names]  
    if("sample_file_regex" %in% names(comp_options)){
      all_files = unique(meta$orig.ident)
      sample_files=all_files[grepl(comp_options$sample_file_regex, all_files)]
      sample_cells = sample_cells[meta[sample_cells, "orig.ident"] %in% sample_files]
    }
    if("sample_files" %in% names(comp_options)){
      sample_files=comp_options$sample_files
      sample_cells = sample_cells[meta[sample_cells, "orig.ident"] %in% sample_files]
    }
    
    control_files=unique(meta[control_cells, "orig.ident"])
    sample_files=unique(meta[sample_cells, "orig.ident"])
    
    if(!all(control_files %in% sample_files)){
      sampleInGroup = 0
    }else{
      sampleInGroup = 1
    }
    
    all_cells<-c(control_cells, sample_cells)
    if("samples" %in% names(comp_options)){
      samples<-comp_options$samples
      all_cells = all_cells[meta[all_cells, "orig.ident"] %in% samples]
    }
    
    de_obj<-subset(obj, cells=all_cells)
    de_obj$Group<-c(rep("control", length(control_cells)), rep("sample", length(sample_cells)))
    de_obj$DisplayGroup<-c(rep(controlGroup, length(control_cells)), rep(sampleGroup, length(sample_cells)))
    
    designdata<-data.frame("Group"=de_obj$Group, "Cell"=colnames(de_obj), "Sample"=de_obj$orig.ident, "DisplayGroup"=de_obj$DisplayGroup)
    designfile<-paste0(prefix, ".design")
    write.csv(designdata, file=designfile, row.names=F, quote=F)
    
    #predefined genes
    if("genes" %in% comp_groups$Key){
      genes_list<-comp_groups$Value[comp_groups$Key=="genes"]
      genes=paste(genes_list, collapse = ",")
    }else{
      genes=""
    }
    
    curdf<-data.frame(prefix=prefix, cellType="", comparison=comp, sampleInGroup=sampleInGroup, design=designfile, genes=genes, stringsAsFactors = F)
    if (is.null(designMatrix)){
      designMatrix = curdf
    }else{
      designMatrix = rbind(designMatrix, curdf)
    }
  }else{
    cts = unique(clusterDf[order(clusterDf$seurat_clusters, decreasing = T), cluster_name])
    prefixList<-gsub(" ", "_", cts)
    prefixList<-gsub(":", "_", prefixList)
    prefixList<-gsub("_+", "_", prefixList)
    
    idx<-1
    for (idx in c(1:length(cts))){
      ct = cts[idx]
      prefix = paste0(prefixList[idx], ".", comp)
      
      clusterCt<-clusterDf[clusterDf[,cluster_name] == ct,]
      de_obj<-subset(obj, cells=rownames(clusterCt))
      clusterCt$sample=de_obj$orig.ident

      invalid_control_names= control_names[!(control_names %in% unique(clusterCt$sample))]
      invalid_sample_names= sample_names[!(sample_names %in% unique(clusterCt$sample))]

      if (length(invalid_control_names) == length(control_names)){
        warning(paste0("There were no control ", paste0(invalid_control_names, collapse=","), " found in cluster ", ct))
        next
      }
      
      if (length(invalid_sample_names)  == length(sample_names)){
        warning(paste0("There were no sample ", paste0(invalid_sample_names, collapse=","), " found in cluster ", ct))
        next
      }
      
      control_cells<-rownames(clusterCt)[clusterCt$sample %in% control_names]  
      sample_cells<-rownames(clusterCt)[clusterCt$sample %in% sample_names]  
      
      all_cells<-c(control_cells, sample_cells)
      
      de_obj<-subset(obj, cells=all_cells)
      de_obj$Group<-c(rep("control", length(control_cells)), rep("sample", length(sample_cells)))
      de_obj$DisplayGroup<-c(rep(controlGroup, length(control_cells)), rep(sampleGroup, length(sample_cells)))
      
      designdata<-data.frame("Group"=de_obj$Group, "Cell"=colnames(de_obj), "Sample"=de_obj$orig.ident, "DisplayGroup"=de_obj$DisplayGroup)
      designfile<-paste0(prefix, ".design")
      write.csv(designdata, file=designfile, row.names=F, quote=F)
      
      curdf<-data.frame(prefix=prefix, cellType=ct, comparison=comp, sampleInGroup=0, design=designfile, stringsAsFactors = F)
      if (is.null(designMatrix)){
        designMatrix = curdf
      }else{
        designMatrix = rbind(designMatrix, curdf)
      }
    }
  }
}

result<-NULL
idx<-1
for(idx in c(1:nrow(designMatrix))){
  prefix=designMatrix[idx, "prefix"]
  designfile=designMatrix[idx, "design"]
  cellType=designMatrix[idx, "cellType"]
  comp=designMatrix[idx, "comparison"]
  sampleInGroup=designMatrix[idx, "sampleInGroup"]
  genes=designMatrix[idx, "genes"]
  if(is.null(genes)){
    genes=""
  }
  
  designdata<-read.csv(designfile, stringsAsFactors = F)
  groups<-designdata$Group
  
  de_obj<-subset(obj, cells=designdata$Cell)
  cells<-as.matrix(de_obj[["RNA"]]@counts)
    
  #filter genes with zero count
  cells<-cells[rowSums(cells)>0,]
  
  #filter genes by tpm
  tpm = sweep(cells, 2, colSums(cells)/1e6, "/")
  min_sample<-filter_cellPercentage * ncol(cells)
  keep_rows <- rowSums(tpm > filter_minTPM) >= min_sample
  rm(tpm)
  
  if(genes != ""){
    gene_list=unlist(strsplit(genes, ','))
    keep2<-rownames(cells) %in% gene_list
    keep_rows = keep_rows | keep2
  }
  
  cells<-cells[keep_rows,]
  cdr <- scale(colMeans(cells > 0))
  
  if(sampleInGroup){
    samples<-designdata$Sample
    design <- model.matrix(~ cdr + samples + groups)
  }else{
    design <- model.matrix(~ cdr + groups)
  }
  
  rownames(design)<-colnames(cells)
  write.csv(design, file=paste0(prefix, ".design_matrix.csv"), quote=F)
  
  cat(prefix, "\n")
  
  dge<-DGEList(cells, group=groups)
  cat("  calcNormFactors", "\n")
  dge<-calcNormFactors(dge)

  cat("  estimateDisp", "\n")
  dge<-estimateDisp(dge,design=design)
  
  if(genes != ""){
    dge<-dge[rownames(dge) %in% gene_list,]
  }
  
  cat("  glmQLFit", "\n")
  fitqlf<-glmQLFit(dge,design=design,robust=TRUE)
  qlf<-glmQLFTest(fitqlf)
  out<-topTags(qlf, n=Inf)
  dge_filename <-paste0(prefix, ".csv")
  write.csv(out$table, file=dge_filename, quote=F)

  if(useRawPvalue){
    sigout<-out$table[(out$table$PValue<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }else{
    sigout<-out$table[(out$table$FDR<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }
  sigFile<-paste0(prefix, ".sig.csv")
  write.csv(sigout, file=sigFile, quote=F)
  
  siggenes<-data.frame(gene=rownames(sigout), stringsAsFactors = F)
  sigGenenameFile<-paste0(prefix, ".sig_genename.txt")
  write.table(siggenes, file=sigGenenameFile, row.names=F, col.names=F, sep="\t", quote=F)
  
  gseaFile<-paste0(prefix, "_GSEA.rnk")
  rankout<-data.frame(gene=rownames(out), sigfvalue=sign(out$table$logFC) * out$table$F)
  write.table(rankout, file=gseaFile, row.names=F, col.names=F, sep="\t", quote=F)
  
  curDF<-data.frame("prefix"=prefix, "cellType"=cellType, "comparison"=comp, "betweenCluster"=bBetweenCluster, "sampleInGroup"=sampleInGroup, "deFile"=dge_filename, "sigFile"=sigFile, "sigGenenameFile"=sigGenenameFile, "gseaFile"=gseaFile, "designFile"=designfile)
  if(is.null(result)){
    result<-curDF
  }else{
    result<-rbind(result, curDF)
  }
}

write.csv(result, file=paste0(outFile, ".edgeR.files.csv"), quote=F)

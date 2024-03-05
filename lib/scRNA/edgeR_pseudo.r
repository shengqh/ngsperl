rm(list=ls()) 
outFile='SADIE_adipose'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/data/shah_lab/shengq2/20240226_mono_scRNA_SADIE/adipose.DE.rds'
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20240304_mona_scRNA_SADIE/20240304_DE_fold1.5/files_edgeR_inCluster_bySample/result')

### Parameter setting end ###

source("scRNA_func.r")
library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(testit)
library(data.table)

options_table<-read.table(parSampleFile3, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
bBetweenCluster<-is_one(myoptions$bBetweenCluster)
pvalue<-as.numeric(myoptions$pvalue)
foldChange<-as.numeric(myoptions$foldChange)
useRawPvalue<-is_one(myoptions$useRawPvalue)
cluster_name=myoptions$cluster_name
min_cell_per_sample=as.numeric(myoptions$filter_min_cell_per_sample)

if(!exists('obj')){
  obj<-read_object(parFile1, parFile2, cluster_name)
  if(!cluster_name %in% colnames(obj@meta.data)){
    if(cluster_name == "bulk"){
      obj=AddMetaData(obj, "bulk", col.name="bulk")
    }
  }
  if(!cluster_name %in% colnames(obj@meta.data)){
    stop(paste0("cluster_name ", cluster_name, " not found in meta.data"))
  }
  obj@meta.data[,cluster_name]<-gsub("^\\s+", "", obj@meta.data[,cluster_name])
  if(!is.null(myoptions$sample_column)){
    if(myoptions$sample_column != ""){
      if(!myoptions$sample_column %in% colnames(obj@meta.data)){
        stop(paste0("sample_column ", myoptions$sample_column, " not found in meta.data"))
      }
      
      obj$orig.ident = obj@meta.data[,myoptions$sample_column]
    }
  }
}

if(myoptions$reduction != "umap"){
  obj[["umap"]] = obj[[myoptions$reduction]]
}

if(1){
  meta<-obj@meta.data
  mt<-data.frame(table(meta[,cluster_name], meta$orig.ident))
  colnames(mt)<-c("cell_type", "sample", "num_cell")
  write.csv(mt, paste0(outFile, ".num_cell.csv"), row.names=F)

  sample_count_df<-get_seurat_sum_count(obj = obj, 
                                        cluster_name = cluster_name, 
                                        min_cell_per_sample = min_cell_per_sample)

  cts<-sample_count_df$cluster
  cts_name_map<-unlist(split(sample_count_df$prefix, sample_count_df$cluster))
  cts_file_map<-unlist(split(sample_count_df$pusedo_file, sample_count_df$cluster))

  comparisons<-read.table(parSampleFile2, stringsAsFactors = F)
  if(ncol(comparisons) == 3){
    colnames(comparisons)<-c("Value", "Key", "Comparison")
  }else{
    colnames(comparisons)<-c("Value", "Comparison")
    comparisons$Key = "groups"
  }

  comparisonNames<-unique(comparisons$Comparison)

  designMatrix<-NULL

  comp <-comparisonNames[1]
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
    
    if("covariances" %in% names(comp_options)){
      if(file.exists(myoptions$covariance_file)){
        covariances_tbl<-data.frame(fread(myoptions$covariance_file, stringsAsFactors = F, header=T), row.names=1)
        assert(all(comp_options$covariances %in% colnames(covariances_tbl)))
        covariances=comp_options$covariances
      }else{
        assert(all(comp_options$covariances %in% colnames(obj@meta.data)))
        covariances_tbl<-unique(obj@meta.data[,c("orig.ident", comp_options$covariances)])
        rownames(covariances_tbl)<-covariances_tbl$orig.ident
        covariances=comp_options$covariances
      }
    }else{
      covariances=NULL
    }

    idx<-1
    for (idx in c(1:length(cts))){
      ct<-cts[idx]
      prefix<-paste0(cts_name_map[ct], ".", comp)
      cts_file<-cts_file_map[ct]
      
      counts<-read.csv(cts_file, row.names = 1)
      cts_control_names<-intersect(control_names, colnames(counts))
      cts_sample_names<-intersect(sample_names, colnames(counts))
      
      if (length(cts_control_names) < 2){
        writeLines(paste0("There were no enough controls: ", paste0(cts_control_names, collapse=",")), paste0(prefix, ".error"))
        next
      }
      
      if (length(cts_sample_names) < 2){
        writeLines(paste0("There were no enough samples: ", paste0(cts_sample_names, collapse=",")), paste0(prefix, ".error"))
        next
      }
      
      design_data<-data.frame("Group"=c(rep("control", length(cts_control_names)), rep("sample", length(cts_sample_names))),
                              "Sample"=c(cts_control_names,cts_sample_names),
                              "DisplayGroup"=c(rep(controlGroup, length(cts_control_names)), rep(sampleGroup, length(cts_sample_names))))
      if(!is.null(covariances)){
        for(cov_name in covariances){
          design_data[,cov_name] = unlist(covariances_tbl[design_data$Sample, cov_name])
        }
      }

      design_file<-paste0(prefix, ".design")
      write.csv(design_data, file=design_file, row.names=F, quote=F)
      
      curdf<-data.frame(prefix=prefix, cellType=ct, comparison=comp, count_file=cts_file, design=design_file, stringsAsFactors = F)
      if (is.null(designMatrix)){
        designMatrix = curdf
      }else{
        designMatrix = rbind(designMatrix, curdf)
      }
    }
  }

  saveRDS(designMatrix, paste0(outFile, ".designMatrix.rds"))
}else{
  designMatrix=readRDS(paste0(outFile, ".designMatrix.rds"))
}

meta<-obj@meta.data

result<-NULL
idx<-18
for(idx in c(1:nrow(designMatrix))){
  prefix=designMatrix[idx, "prefix"]
  design_file=designMatrix[idx, "design"]
  cellType=designMatrix[idx, "cellType"]
  comp=designMatrix[idx, "comparison"]
  count_file=designMatrix[idx,"count_file"]

  design_data<-read.csv(design_file, stringsAsFactors = F)

  covariances=colnames(design_data)[!(colnames(design_data) %in% c("Group", "Sample", "DisplayGroup"))]

  while(TRUE){
    if(length(covariances) > 0){
      variables = c(covariances, "Group")
    }else{
      variables = c("Group")
    }
    formula_str = paste0("~ ", paste0(variables, collapse = " + "))
    design <- model.matrix(as.formula(formula_str), design_data)

    rownames(design)<-design_data$Sample
    write.csv(design, file=paste0(prefix, ".design_matrix.csv"), quote=F)
    
    cat(prefix, "\n")

    counts<-read.csv(count_file, row.names = 1)
    counts<-counts[,design_data$Sample]
    counts_file = paste0(prefix, ".counts.csv")
    write.csv(counts, file=counts_file, quote=F)

    groups<-design_data$Group
    
    dge <- DGEList(counts=counts, group=groups)
    keep <- filterByExpr(dge, group=groups)
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    
    cat("  calcNormFactors", "\n")
    dge<-calcNormFactors(dge,method = "TMM")

    cat("  estimateDisp", "\n")
    dge<-estimateDisp(dge,design=design)

    if("common.dispersion" %in% names(dge)){
      if(!is.na(dge$common.dispersion)){
        break
      }
    }

    if(length(covariances) == 0){
      stop("estimateDisp failed")
    }

    covariances = c()
    cat("  redo estimateDisp without covariates\n")
    design_data = design_data[,colnames(design_data) %in% c("Group", "Sample", "DisplayGroup")]
    write.csv(design_data, design_file, row.names=FALSE)
  }
  
  saveRDS(dge, paste0(prefix, ".dge.rds"))

  #https://www.nature.com/articles/s41467-021-25960-2#Sec39
  cat("  glmFit", "\n")
  fit<-glmFit(dge,design=design,robust=TRUE)
  fitTest<-glmLRT(fit)
  out<-topTags(fitTest, n=Inf)
  dge_filename <-paste0(prefix, ".csv")
  write.csv(out$table, file=dge_filename, quote=F)

  log_cpm <- cpm(dge, log=TRUE)
  cpm_file = paste0(prefix, ".cpm.csv")
  write.csv(log_cpm, file=cpm_file, quote=F)

  if(useRawPvalue){
    sigout<-out$table[(out$table$PValue<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }else{
    sigout<-out$table[(out$table$FDR<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }
  sigFile<-paste0(prefix, ".sig.csv")
  write.csv(sigout, file=sigFile, quote=F)
  
  sig_genes<-data.frame(gene=rownames(sigout), stringsAsFactors = F)
  sigGenenameFile<-paste0(prefix, ".sig_genename.txt")
  write.table(sig_genes, file=sigGenenameFile, row.names=F, col.names=F, sep="\t", quote=F)
  
  gseaFile<-paste0(prefix, "_GSEA.rnk")
  rankout<-data.frame(gene=rownames(out), sigfvalue=sign(out$table$logFC) * (-log10(out$table$PValue)))
  write.table(rankout, file=gseaFile, row.names=F, col.names=F, sep="\t", quote=F)

  if(nrow(sig_genes) > 0){
    sig_gene = sig_genes$gene[1]

    ct_meta<-meta[meta[,cluster_name] == cellType,]
    ct_cells<-rownames(ct_meta)[ct_meta$orig.ident %in% design_data$Sample]
    cell_obj<-subset(obj, cells=ct_cells)

    g<-get_sig_gene_figure(cell_obj, sigout, design_data, sig_gene, DE_by_cell=FALSE, is_between_cluster=FALSE, log_cpm=log_cpm)
    png(paste0(prefix, ".top_1_gene.png"), width=3000, height=2500, res=300)
    print(g)
    dev.off()
  }
  
  curDF<-data.frame("prefix"=prefix, "cellType"=cellType, "comparison"=comp, "betweenCluster"=0, "sampleInGroup"=0, "deFile"=dge_filename, "sigFile"=sigFile, "sigGenenameFile"=sigGenenameFile, "gseaFile"=gseaFile, "designFile"=design_file, "cpmFile"=cpm_file)
  if(is.null(result)){
    result<-curDF
  }else{
    result<-rbind(result, curDF)
  }
}

write.csv(result, file=paste0(outFile, ".edgeR.files.csv"), quote=F)

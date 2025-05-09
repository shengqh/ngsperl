rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose/result/combined.final.rds'
parFile2='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose/result/combined.meta.rds'
parFile3=''


setwd('/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose_edgeR_inCluster_byCell/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
source("scRNA_func.r")
library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(testit)

MIN_NUM_CELL=10

sumcount<-function(ct_count, names, sample_df){
  result<-lapply(names, function(x){
    res<-ct_count[,sample_df$Cell[sample_df$Sample ==x],drop=F]
    apply(res, 1, sum)
  })
  rescount<-do.call(cbind, result)
  colnames(rescount)<-names
  return(rescount)
}

options_table<-read.table(parSampleFile3, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
bBetweenCluster<-ifelse(myoptions$bBetweenCluster == "0", FALSE, TRUE)
filter_cellPercentage<-as.numeric(myoptions$filter_cellPercentage)
filter_minTPM<-as.numeric(myoptions$filter_minTPM)
pvalue<-as.numeric(myoptions$pvalue)
foldChange<-as.numeric(myoptions$foldChange)
useRawPvalue<-ifelse(myoptions$useRawPvalue == "0", FALSE, TRUE)
cluster_name=myoptions$cluster_name

if(!exists('obj')){
  obj<-read_object(parFile1, parFile2, cluster_name)
}

detail_folder = paste0(outFile, ".edgeR_by_cell/")
if(!dir.exists(detail_folder)){
  dir.create(detail_folder)
}
detail_prefix = paste0(detail_folder, outFile)

meta<-obj@meta.data
# if("seurat_cell_type" %in% colnames(meta)){
#   mt<-data.frame(table(meta$seurat_cell_type, meta$orig.ident))
# }else{
mt<-data.frame(table(meta[,cluster_name], meta$orig.ident))
# }
colnames(mt)<-c("cell_type", "sample", "num_cell")
write.csv(mt, paste0(detail_prefix, ".num_cell.csv"), row.names=F)

clusterDf<-obj@meta.data

comparisons<-read.table(parSampleFile2, sep="\t", stringsAsFactors = F, fill=TRUE, header=F)
if(ncol(comparisons) == 3){
  colnames(comparisons)<-c("Value", "Key", "Comparison")
}else{
  colnames(comparisons)<-c("Value", "Comparison")
  comparisons$Key = "groups"
}

comparisonNames<-unique(comparisons$Comparison)

comp <-comparisonNames[1]

cat("start to generate design files\n")
designFailed<-data.frame("comp"=character(), "celltype"=character(), "reason"=character())
designMatrix<-NULL
for (comp in comparisonNames){
  cat(comp, "\n")
  comp_groups<-comparisons[comparisons$Comparison==comp,]
  comp_options = split(comp_groups$Value, comp_groups$Key)
  
  if("groups" %in% names(comp_options)){
    sampleGroups<-read.table(parSampleFile1, sep="\t", stringsAsFactors = F)
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
    covariances_tbl<-read.table(myoptions$covariance_file, sep="\t", stringsAsFactors = F, header=T)
    assert(all(comp_options$covariances %in% colnames(covariances_tbl)))
    for(col_name in comp_options$covariances){
      colmap = unlist(split(covariances_tbl[,col_name], covariances_tbl$Sample))
      obj<-AddMetaData(obj, colmap[obj$orig.ident], col.name = col_name )
    }
    covariances=comp_options$covariances
  }else{
    covariances=NULL
  }

  if(bBetweenCluster){
    prefix<-get_valid_path(comp)
    file_prefix = paste0(detail_prefix, ".", prefix)
    
    control_cells<-rownames(meta)[meta[,cluster_name] %in% control_names]  
    if (length(control_cells) == 0){
      stop(paste0("no control cells found, check your control names:", paste0(control_names, collapse = ",")))
    }
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
    if (length(sample_cells) == 0){
      stop(paste0("no sample cells found, check your sample names:", paste0(sample_cells, collapse = ",")))
    }
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
      sampleInGroup = FALSE
    }else{
      sampleInGroup = TRUE
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
    
    designfile<-paste0(file_prefix, ".design")
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
    if(myoptions$DE_cluster_pattern != '*' & myoptions$DE_cluster_pattern != ''){
      cts=cts[grepl(myoptions$DE_cluster_pattern, cts)]
    }
    
    prefixList<-get_valid_path(cts)

    idx<-5
    for (idx in c(1:length(cts))){
      ct = cts[idx]
      cat("  ", as.character(ct), "\n")
      prefix = paste0(prefixList[idx], ".", comp)
      ct_file_prefix = paste0(detail_prefix, ".", prefix)
      
      clusterCt<-clusterDf[clusterDf[,cluster_name] == ct,]
      de_obj<-subset(obj, cells=rownames(clusterCt))
      clusterCt$sample=de_obj$orig.ident

      invalid_control_names= control_names[!(control_names %in% unique(clusterCt$sample))]
      invalid_sample_names= sample_names[!(sample_names %in% unique(clusterCt$sample))]

      if (length(invalid_control_names) == length(control_names)){
        error_msg = paste0("There were no control ", paste0(invalid_control_names, collapse=","), " found in object sample names!")
        designFailed[nrow(designFailed) + 1,] <- c(comp, as.character(ct), error_msg)
        cat(error_msg, "\n", file=stderr())
        next
      }
      
      if (length(invalid_sample_names)  == length(sample_names)){
        error_msg = paste0("There were no sample ", paste0(invalid_sample_names, collapse=","), " found in object sample names!")
        designFailed[nrow(designFailed) + 1,] <- c(comp, as.character(ct), error_msg)
        cat(error_msg, "\n", file=stderr())
        next
      }
      
      control_cells<-rownames(clusterCt)[clusterCt$sample %in% control_names]  
      sample_cells<-rownames(clusterCt)[clusterCt$sample %in% sample_names]  

      if(length(control_cells) < MIN_NUM_CELL){
        error_msg = paste0("There were only ", length(control_cells), " cells found in control group, less than required ", MIN_NUM_CELL, "!")
        designFailed[nrow(designFailed) + 1,] <- c(comp, as.character(ct), error_msg)
        cat(error_msg, "\n", file=stderr())
        next
      }

      if(length(sample_cells) < MIN_NUM_CELL){
        error_msg = paste0("There were only ", length(sample_cells), " cells found in sample group, less than required ", MIN_NUM_CELL, "!")
        designFailed[nrow(designFailed) + 1,] <- c(comp, as.character(ct), error_msg)
        cat(error_msg, "\n", file=stderr())
        next
      }

      all_cells<-c(control_cells, sample_cells)
      
      de_obj<-subset(obj, cells=all_cells)
      de_obj$Group<-ifelse(colnames(de_obj) %in% control_cells, "control", "sample")
      de_obj$DisplayGroup<-ifelse(colnames(de_obj) %in% control_cells, controlGroup, sampleGroup)
      
      designdata<-data.frame("Group"=de_obj$Group, "Cell"=colnames(de_obj), "Sample"=de_obj$orig.ident, "DisplayGroup"=de_obj$DisplayGroup)
      if(!is.null(covariances)){
        for(cov_name in covariances){
          designdata[,cov_name] = unlist(de_obj@meta.data[,cov_name])
        }
      }
      assert(all(covariances %in% colnames(designdata)))

      designfile<-paste0(ct_file_prefix, ".design")
      write.csv(designdata, file=designfile, row.names=F, quote=F)

      designdatatbl=designdata |>
        dplyr::group_by(Group, DisplayGroup, Sample) |>
        dplyr::summarize(num_cell=n(), .groups = "drop")
      write.csv(designdatatbl, paste0(designfile, ".num_cell.csv"), row.names=F, quote=F)
      
      curdf<-data.frame(prefix=prefix, cellType=ct, comparison=comp, sampleInGroup=FALSE, design=designfile, stringsAsFactors = F)
      if (is.null(designMatrix)){
        designMatrix = curdf
      }else{
        designMatrix = rbind(designMatrix, curdf)
      }
    }
  }
}

if(nrow(designFailed) > 0){
  write.csv(designFailed, paste0(detail_prefix, ".design_failed.csv"), row.names=F)
}

design_matrix_file=paste0(detail_prefix, ".design_matrix.csv")
write.csv(designMatrix, file=design_matrix_file, row.names=F)
designMatrix=read.csv(design_matrix_file, stringsAsFactors = F)

cat("start to perform edgeR\n")
result<-NULL
idx<-1
for(idx in c(1:nrow(designMatrix))){
  prefix=designMatrix[idx, "prefix"]
  file_prefix = paste0(detail_prefix, ".", prefix)

  cat("performing", prefix, "...\n")

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
  #designdata is generated from corresponding de_obj, so the order of cells should be the same
  stopifnot(colnames(de_obj) == designdata$Cell)

  cells<-as.matrix( MyGetAssayData(de_obj, "RNA", "counts"))

  #filter genes with zero count
  cells<-cells[rowSums(cells)>0,]

  control_cells=cells[,designdata[designdata$Group=="control", "Cell"]]
  sample_cells=cells[,designdata[designdata$Group=="sample", "Cell"]]

  do_filter<-function(cur_cells, filter_cellPercentage, filter_minTPM, min_sample){
    #filter genes by tpm
    tpm = sweep(cur_cells, 2, colSums(cur_cells)/1e6, "/")
    min_sample<-filter_cellPercentage * ncol(cur_cells)
    keep_rows <- rowSums(tpm > filter_minTPM) >= min_sample
    return(keep_rows)
  }
  filter_control=do_filter(control_cells, filter_cellPercentage, filter_minTPM, min_sample)
  filter_sample=do_filter(sample_cells, filter_cellPercentage, filter_minTPM, min_sample)
  keep_rows = filter_control | filter_sample
  
  if(genes != ""){
    gene_list=unlist(strsplit(genes, ','))
    keep2<-rownames(cells) %in% gene_list
    keep_rows = keep_rows | keep2
  }
  
  cells<-cells[keep_rows,]
  cdr <- scale(colMeans(cells > 0))
  designdata$cdr = unlist(cdr)
  
  variables = c("cdr")
  if(sampleInGroup){
    # For comparison in each cluster, the sample would be in each group, if we put sample in the design matrix, 
    # it will be confounded with group and failed in estimateDisp: 
    # Design matrix not of full rank.  The following coefficients not estimable: Groupsample
    if(length(unique(designdata$Sample)) > 1){
      variables = c(variables, "Sample")
    }
  }
  if(!is.null(covariances)){
    variables = c(variables, covariances)
  }
  variables = c(variables, "Group")
  formula_str = paste0("~ ", paste0(variables, collapse = " + "))
  cat(formula_str, "\n")

  design <- model.matrix(as.formula(formula_str), designdata)

  rownames(design)<-colnames(cells)
  write.csv(design, file=paste0(file_prefix, ".design_matrix.csv"), quote=F)
  
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

  cat("  glmQLFTest", "\n")
  qlf<-glmQLFTest(fitqlf)

  cat("  topTags", "\n")
  out<-topTags(qlf, n=Inf)
  dge_filename <-paste0(file_prefix, ".csv")
  write.csv(out$table, file=dge_filename, quote=F)

  cat("  volcano plot", "\n")
  save_volcano_plot(edgeR_out_table=out$table, 
                                  prefix=file_prefix, 
                                  useRawPvalue=useRawPvalue, 
                                  pvalue=pvalue, 
                                  foldChange=foldChange, 
                                  comparisonTitle=paste0(cellType, " : ", comp))

  if(useRawPvalue){
    sigout<-out$table[(out$table$PValue<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }else{
    sigout<-out$table[(out$table$FDR<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }
  sigFile<-paste0(file_prefix, ".sig.csv")
  write.csv(sigout, file=sigFile, quote=F)

  if(nrow(sigout) > 0){
    sig_gene=rownames(sigout)[1]
    g<-get_sig_gene_figure(de_obj, sigout, designdata, sig_gene, DE_by_cell=TRUE, is_between_cluster=bBetweenCluster, log_cpm=NULL)
    ggsave(paste0(file_prefix, ".top_1_gene.png"),  g, width=3000, height=2500, units="px", dpi=300)
  }
  
  siggenes<-data.frame(gene=rownames(sigout), stringsAsFactors = F)
  sigGenenameFile<-paste0(file_prefix, ".sig_genename.txt")
  write.table(siggenes, file=sigGenenameFile, row.names=F, col.names=F, sep="\t", quote=F)
  
  gseaFile<-paste0(file_prefix, "_GSEA.rnk")
  pValuesNoZero=out$table$PValue
  if (any(pValuesNoZero == 0)){
    pValuesNoZero[pValuesNoZero == 0] = min(pValuesNoZero[pValuesNoZero > 0],na.rm=TRUE)/10
  }
  rankout<-data.frame(gene=rownames(out), sigfvalue=sign(out$table$logFC) * (-log10(pValuesNoZero)))
  rankout<-rankout[order(rankout$sigfvalue, decreasing=TRUE),]
  write.table(rankout, file=gseaFile, row.names=F, col.names=F, sep="\t", quote=F)
  
  curDF<-data.frame("prefix"=prefix, "cellType"=cellType, "comparison"=comp, "betweenCluster"=bBetweenCluster, "sampleInGroup"=sampleInGroup, "deFile"=dge_filename, "sigFile"=sigFile, "sigGenenameFile"=sigGenenameFile, "gseaFile"=gseaFile, "designFile"=designfile)
  if(is.null(result)){
    result<-curDF
  }else{
    result<-rbind(result, curDF)
  }
  cat("  done", "\n")
}

cat("all done\n")
write.csv(result, file=paste0(outFile, ".edgeR.files.csv"), quote=F)

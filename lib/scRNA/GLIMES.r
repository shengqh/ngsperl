rm(list=ls()) 
outFile='monocytes'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/data/wanjalla_lab/projects/20250420_P12891-P12795_10Flex_hg38_cellbender/20260821_PADI2/20260821_PADI2_01_analysis/monocyte_obj.rds'
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20250420_P12891-P12795_10Flex_hg38_cellbender/20260821_PADI2/20260821_PADI2_02_DE/final_obj_GLIMES_betweenCluster_byCell/result')

### Parameter setting end ###

source("scRNA_func.r")
source("countTableVisFunctions.R")
library("ggplot2")
library("ggpubr")
library("Seurat")
library("testit")
library("GLIMES")
library("EnhancedVolcano")

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
glmm_method=myoptions$glmm_method

if(!exists('obj')){
  obj<-read_object(parFile1, parFile2, cluster_name)
}

detail_folder = paste0(outFile, ".GLIMES_by_cell/")
if(!dir.exists(detail_folder)){
  dir.create(detail_folder)
}
detail_prefix = paste0(detail_folder, outFile)

meta<-obj@meta.data

mt<-data.frame(table(meta[,cluster_name], meta$orig.ident))

colnames(mt)<-c("cluster", "sample", "num_cell")
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
    groups<-comp_options$groups
    controlGroup<-groups[1]
    sampleGroup<-groups[2]

    if(bBetweenCluster){
      control_names<-controlGroup
      sample_names<-sampleGroup
    }else{
      sampleGroups<-read.table(parSampleFile1, sep="\t", stringsAsFactors = F)
      colnames(sampleGroups)<-c("Sample","Group")
      
      control_names<-sampleGroups$Sample[sampleGroups$Group==controlGroup]
      sample_names<-sampleGroups$Sample[sampleGroups$Group==sampleGroup]
    }
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

    sample_cells<-rownames(meta)[meta[,cluster_name] %in% sample_names]  
    if (length(sample_cells) == 0){
      stop(paste0("no sample cells found, check your sample names:", paste0(sample_cells, collapse = ",")))
    }
    
    all_cells<-c(control_cells, sample_cells)
    
    de_obj<-subset(obj, cells=all_cells)
    de_obj@meta.data$Group<-ifelse(colnames(de_obj) %in% control_cells, "control", "sample")
    de_obj@meta.data$DisplayGroup<-ifelse(colnames(de_obj) %in% control_cells, controlGroup, sampleGroup)
    
    designdata<-data.frame( "Group"=de_obj$Group, 
                            "Cell"=colnames(de_obj), 
                            "Sample"=de_obj$orig.ident, 
                            "DisplayGroup"=de_obj$DisplayGroup)

    designfile<-paste0(file_prefix, ".design")
    write.csv(designdata, file=designfile, row.names=F, quote=F)
    
    #predefined genes
    if("genes" %in% comp_groups$Key){
      genes_list<-comp_groups$Value[comp_groups$Key=="genes"]
      genes=paste(genes_list, collapse = ",")
    }else{
      genes=""
    }
    
    curdf<-data.frame(prefix=prefix, cellType="", comparison=comp, sampleInGroup=FALSE, design=designfile, genes=genes, stringsAsFactors = F)
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
      de_obj@meta.data$Group<-ifelse(colnames(de_obj) %in% control_cells, "control", "sample")
      de_obj@meta.data$DisplayGroup<-ifelse(colnames(de_obj) %in% control_cells, controlGroup, sampleGroup)
      
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

cat("start to perform GLIMES\n")
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
  if(is.na(genes)){
    genes=""
  }

  if(is.null(cellType)){
    cellType=""
  }
  if(is.na(cellType)){
    cellType=""
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
  #remove cells with too few total reads after filtering and also change designdata and group
  cellTotalReads=colSums(cells)
  if (any(cellTotalReads<50)) {
    selectedCellsInd=which(cellTotalReads>=50)
    cells=cells[,selectedCellsInd]
    groups=groups[selectedCellsInd]
    designdata=designdata[selectedCellsInd,]
  }

  variables=c()
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
  
  de_obj = subset(de_obj, cells=designdata$Cell)
  stopifnot(colnames(de_obj) == designdata$Cell)

  de_obj@meta.data$Group = factor(designdata$Group, levels=c("control", "sample"))
  de_obj@meta.data$Sample = designdata$Sample
  de_obj@meta.data$DisplayGroup = designdata$DisplayGroup

  de_sce = as.SingleCellExperiment(de_obj)

  uniq_groups = unique(de_obj@meta.data |> dplyr::select(Group, Sample))

  if(any(table(uniq_groups$Group) == 1)){ # if there is group with only one sample, replicates should be all same.
    replicates = ct_col
  }else{
    replicates = "Sample"
  }
  
  if(glmm_method == "binomial"){
    cat("Using binomial GLMM for DE analysis\n")
    glmm_df = binomial_glmm_DE(sce = de_sce, comparison = "Group", replicates = replicates, other_fixed=covariances)
  }else{
    cat("Using poisson GLMM for DE analysis\n")
    glmm_df = poisson_glmm_DE(sce = de_sce, comparison = "Group", replicates = "Sample", other_fixed=covariances)
  }

  if(useRawPvalue){
    glmm_df$DEGs = identifyDEGs(adj_pval = glmm_df$pval, 
                                log2FC = glmm_df$log2FC, 
                                log2mean = glmm_df$log2mean, 
                                log2meandiff = glmm_df$log2meandiff, 
                                log2FCcutoff = log2(1.2),
                                newcriteria = T)
  }else{
    glmm_df$DEGs = identifyDEGs(adj_pval = glmm_df$BH, 
                                log2FC = glmm_df$log2FC, 
                                log2mean = glmm_df$log2mean, 
                                log2meandiff = glmm_df$log2meandiff, 
                                log2FCcutoff = log2(1.2),
                                newcriteria = T)
  }
  glmm_df = glmm_df[order(glmm_df$pval, decreasing = F),]

  #change to edgeR compatible column names
  glmm_df=glmm_df |>
    dplyr::mutate(PValue=pval, FDR=BH, logFC=log2FC) 

  dge_filename <-paste0(file_prefix, ".csv")
  write.csv(glmm_df, file=dge_filename, quote=F, row.names=FALSE)
  #glmm_df=read.csv(dge_filename, stringsAsFactors = F)

  glmm_df <- glmm_df |> dplyr::filter(!is.na(PValue))

  sigout = glmm_df |> dplyr::filter(DEGs == TRUE)
  sigFile<-paste0(file_prefix, ".sig.csv")
  write.csv(sigout, file=sigFile, quote=F)

  # filter extreme logFC values to avoid plotting issues
  if(min(glmm_df$logFC) < -20){
    glmm_df$logFC[glmm_df$logFC < -20] = min(glmm_df$logFC[glmm_df$logFC > -20])
  }
  if(min(glmm_df$logFC) > 20){
    glmm_df$logFC[glmm_df$logFC > 20] = max(glmm_df$logFC[glmm_df$logFC < 20])
  }

  cat("  volcano plot", "\n")
  comp_str=gsub("_vs_", " vs ", comp)

  glmm_df$Status <- ifelse(glmm_df$DEGs, "Significant", "NS")
  keyvals <- ifelse(
    glmm_df$DEGs,
    "red",
    "grey70"
  )
  names(keyvals) <- glmm_df$Status
  genes=glmm_df$genes
  genes[!glmm_df$DEGs] = ""

  yname=bquote(-log[10](p~value))

  g = EnhancedVolcano(
    glmm_df,
    lab = genes,
    x = "logFC",
    y = "PValue",
    title = ifelse(cellType=="", comp_str, paste0(cellType, " : ", comp_str)),
    colCustom = keyvals,
    pCutoff = 1,      # disable internal significance filtering
    FCcutoff = 0,     # disable internal FC filtering
    pointSize = 1,
    labSize = 3.0,
    colAlpha = 1,
    subtitle = NULL
  ) + ylab(yname) + theme(plot.title = element_text(hjust = 0.5))

  for(extension in c("png", "pdf")){
    volcano_file=paste0(file_prefix, ".volcano.", extension)
    ggsave(volcano_file, g, width=7, height=7, units="in", dpi=300, bg="white")
  }
  
  # save_volcano_plot(edgeR_out_table=glmm_df, 
  #                   prefix=file_prefix, 
  #                   useRawPvalue=useRawPvalue, 
  #                   pvalue=pvalue, 
  #                   foldChange=foldChange, 
  #                   comparisonTitle=ifelse(cellType=="", comp_str, paste0(cellType, " : ", comp_str)),
  #                   genes=glmm_df$genes)

  sigout = glmm_df |> dplyr::filter(DEGs == TRUE)

  if(nrow(sigout) > 0) {
    save_top_gene_figures(de_obj, sigout, designdata, bBetweenCluster, file_prefix)
  }
  
  siggenes<-data.frame(gene=sigout$genes, stringsAsFactors = F)
  sigGenenameFile<-paste0(file_prefix, ".sig_genename.txt")
  write.table(siggenes, file=sigGenenameFile, row.names=F, col.names=F, sep="\t", quote=F)
  
  gseaFile<-paste0(file_prefix, "_GSEA.rnk")
  rankout<-data.frame(gene=glmm_df$genes, sigfvalue=sign(glmm_df$logFC) * (-log10(glmm_df$PValue)))
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
write.csv(result, file=paste0(outFile, ".GLIMES.files.csv"), quote=F)

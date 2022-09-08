#rm(list=ls()) 
outFile='AK6383'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose/result/AK6383.final.rds'
parFile2='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose/result/AK6383.meta.rds'
parFile3=''


setwd('C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose_edgeR_inCluster_bySample/result')

### Parameter setting end ###

source("scRNA_func.r")
library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(testit)

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
sample_count_df<-get_seurat_sum_count(obj, cluster_name)

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
  
  if("covariances" %in% names(comp_options)){
    covariances_tbl<-read.table(myoptions$covariance_file, sep="\t", stringsAsFactors = F, header=T, row.names=1)
    assert(all(comp_options$covariances %in% colnames(covariances_tbl)))
    covariances=comp_options$covariances
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
    
    designdata<-data.frame("Group"=c(rep("control", length(cts_control_names)), rep("sample", length(cts_sample_names))),
                           "Sample"=c(cts_control_names,cts_sample_names),
                           "DisplayGroup"=c(rep(controlGroup, length(cts_control_names)), rep(sampleGroup, length(cts_sample_names))))
    if(!is.null(covariances)){
      for(cov_name in covariances){
        designdata[,cov_name] = unlist(covariances_tbl[designdata$Sample, cov_name])
      }
    }

    designfile<-paste0(prefix, ".design")
    write.csv(designdata, file=designfile, row.names=F, quote=F)
    
    curdf<-data.frame(prefix=prefix, cellType=ct, comparison=comp, count_file=cts_file, design=designfile, stringsAsFactors = F)
    if (is.null(designMatrix)){
      designMatrix = curdf
    }else{
      designMatrix = rbind(designMatrix, curdf)
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
  count_file=designMatrix[idx,"count_file"]

  designdata<-read.csv(designfile, stringsAsFactors = F)

  if(!is.null(covariances)){
    variables = c(covariances, "Group")
  }else{
    variables = c("Group")
  }
  formula_str = paste0("~ ", paste0(variables, collapse = " + "))
  design <- model.matrix(as.formula(formula_str), designdata)

  rownames(design)<-designdata$Sample
  write.csv(design, file=paste0(prefix, ".design_matrix.csv"), quote=F)
  
  cat(prefix, "\n")

  counts<-read.csv(count_file, row.names = 1)
  counts<-counts[,designdata$Sample]
  counts_file = paste0(prefix, ".counts.csv")
  write.csv(counts, file=counts_file, quote=F)

  groups<-designdata$Group
  
  dge <- DGEList(counts=counts, group=groups)
  keep <- filterByExpr(dge, group=groups)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  
  cat("  calcNormFactors", "\n")
  dge<-calcNormFactors(dge,method = "TMM")

  cat("  estimateDisp", "\n")
  dge<-estimateDisp(dge,design=design)
  
  saveRDS(dge, paste0(prefix, ".dge.rds"))

  #https://www.nature.com/articles/s41467-021-25960-2#Sec39
  cat("  glmFit", "\n")
  fit<-glmFit(dge,design=design,robust=TRUE)
  fitTest<-glmLRT(fit)
  out<-topTags(fitTest, n=Inf)
  dge_filename <-paste0(prefix, ".csv")
  write.csv(out$table, file=dge_filename, quote=F)

  logcpm <- cpm(dge, log=TRUE)
  cpm_file = paste0(prefix, ".cpm.csv")
  write.csv(logcpm, file=cpm_file, quote=F)

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
  
  curDF<-data.frame("prefix"=prefix, "cellType"=cellType, "comparison"=comp, "betweenCluster"=0, "sampleInGroup"=0, "deFile"=dge_filename, "sigFile"=sigFile, "sigGenenameFile"=sigGenenameFile, "gseaFile"=gseaFile, "designFile"=designfile, "cpmFile"=cpm_file)
  if(is.null(result)){
    result<-curDF
  }else{
    result<-rbind(result, curDF)
  }
}

write.csv(result, file=paste0(outFile, ".edgeR.files.csv"), quote=F)

rm(list=ls()) 
outFile='carotid'
parSampleFile1=''
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/nobackup/shah_lab/shengq2/20240208_CAC_proteomics_scRNA/chiara_scRNA/T01_prepare_data/carotid.DE.rds'
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20240208_CAC_proteomics_scRNA/chiara_scRNA/T03_DE_fold1.2_carotid.no_CV7209/bulk_edgeR_inCluster_bySample/result')

### Parameter setting end ###

source("scRNA_func.r")
source("countTableVisFunctions.R")

library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(testit)
library(data.table)
library(matrixStats)
library(ggrepel)

options_table<-read.table(parSampleFile3, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
bBetweenCluster<-is_one(myoptions$bBetweenCluster)
pvalue<-as.numeric(myoptions$pvalue)
foldChange<-as.numeric(myoptions$foldChange)
useRawPvalue<-is_one(myoptions$useRawPvalue)
cluster_name=myoptions$cluster_name
min_cell_per_sample=as.numeric(myoptions$filter_min_cell_per_sample)

group_column=myoptions$group_column

discard_samples=unlist(strsplit(myoptions$discard_samples, ','))

comparisons<-read.table(parSampleFile2, sep="\t", stringsAsFactors = F, fill=TRUE, header=F)
if(ncol(comparisons) == 3){
  colnames(comparisons)<-c("Value", "Key", "Comparison")
}else{
  colnames(comparisons)<-c("Value", "Comparison")
  comparisons$Key = "groups"
}

comparisonNames<-unique(comparisons$Comparison)

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

  if(length(discard_samples) > 0){
    cat("discard_samples: ", paste0(discard_samples, collapse=", "), "\n")
    discard_cells = colnames(obj)[obj@meta.data$orig.ident %in% discard_samples]
    obj<-subset(obj, cells=discard_cells, invert=TRUE)
  }
}

if(myoptions$reduction != "umap"){
  obj[["umap"]] = obj[[myoptions$reduction]]
}

detail_folder = paste0(outFile, ".edgeR_by_sample/")
if(!dir.exists(detail_folder)){
  dir.create(detail_folder)
}

detail_prefix = paste0(detail_folder, outFile)

if(1){
  meta<-obj@meta.data
  mt<-data.frame(table(meta[,cluster_name], meta$orig.ident))
  colnames(mt)<-c("cell_type", "sample", "num_cell")
  write.csv(mt, paste0(detail_prefix, ".num_cell.csv"), row.names=F)

  sample_count_df<-get_seurat_sum_count(obj = obj, 
                                        cluster_name = cluster_name, 
                                        min_cell_per_sample = min_cell_per_sample,
                                        target_folder = detail_folder)

  cts<-sample_count_df$cluster
  cts_name_map<-unlist(split(sample_count_df$prefix, sample_count_df$cluster))
  cts_file_map<-unlist(split(sample_count_df$pseudo_file, sample_count_df$cluster))

  designMatrix<-NULL

  comp <-comparisonNames[1]
  for (comp in comparisonNames){
    comp_groups<-comparisons[comparisons$Comparison==comp,]
    comp_options = split(comp_groups$Value, comp_groups$Key)
    
    if("groups" %in% names(comp_options)){
      if(!is.null(myoptions$group_column)){
        if(myoptions$group_column != ""){
          if(!myoptions$group_column %in% colnames(obj@meta.data)){
            stop(paste0("group_column ", myoptions$group_column, " not found in meta.data of ", parFile1))
          }
        }

        sampleGroups = unique(obj@meta.data[,c("orig.ident", myoptions$group_column)]) %>%
          dplyr::rename(Sample=orig.ident, Group=!!myoptions$group_column) %>%
          tibble::remove_rownames() %>%
          dplyr::arrange(Group, Sample)
        write.table(sampleGroups, file="fileList1.txt", sep="\t", row.names=F, col.names=F, quote=F)
      }else{
        sampleGroups<-read.table(parSampleFile1, stringsAsFactors = F)
        colnames(sampleGroups)<-c("Sample","Group")
      }
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
      prefix<-paste0(detail_folder, cts_name_map[ct], ".", comp)
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

  saveRDS(designMatrix, paste0(detail_prefix, ".designMatrix.rds"))
}else{
  designMatrix=readRDS(paste0(detail_prefix, ".designMatrix.rds"))
}

drawHCA<-function(prefix, top_mat, design_data, group_colors, usePearsonInHCA=FALSE){
  mat_scaled = t(scale(t(top_mat)))

  ha=HeatmapAnnotation( Group=design_data$DisplayGroup,
                        col=list(Group=group_colors),
                        annotation_legend_param = list(Group = list(ncol = 1, title = "Group", title_position = "topleft")))

  if(usePearsonInHCA){
    clustering_distance_columns="pearson"
  }else{
    clustering_distance_columns="euclidean"
  }

  draw_heatmap_png( filepath=paste0(prefix, ".heatmap.png"), 
                    htdata=mat_scaled, 
                    name="zscore", 
                    show_row_names=FALSE, 
                    show_column_names=TRUE,
                    save_rds=FALSE,
                    clustering_distance_columns=clustering_distance_columns,
                    top_annotation=ha,
                    show_row_dend=FALSE)
}

meta<-obj@meta.data

result<-NULL
idx<-1
for(idx in c(1:nrow(designMatrix))){
  prefix=designMatrix[idx, "prefix"]
  design_file=designMatrix[idx, "design"]
  cellType=designMatrix[idx, "cellType"]
  comp=designMatrix[idx, "comparison"]
  count_file=designMatrix[idx,"count_file"]

  design_data<-read.csv(design_file, stringsAsFactors = F)

  covariances=colnames(design_data)[!(colnames(design_data) %in% c("Group", "Sample", "DisplayGroup"))]

  design_groups=unique(design_data[,c("Group", "DisplayGroup")])
  group_colors=c("blue", "red")
  names(group_colors)<-design_groups$DisplayGroup
  sample_colors=group_colors[design_data$DisplayGroup]
  names(sample_colors)=design_data$Sample

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
    
    dge_all <- DGEList(counts=counts, group=groups)

    #for DE analysis, we need to filter the genes with low counts
    keep <- filterByExpr(dge_all, group=groups)
    dge <- dge_all[keep, , keep.lib.sizes=FALSE]
    
    cat("  calcNormFactors", "\n")
    dge<-calcNormFactors(dge,method = "TMM")

    cat("  plotMDS", "\n")
    mds <- plotMDS(dge, plot=FALSE)
    toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Sample = design_data$Sample, Group = design_data$DisplayGroup)
    g<-ggplot(toplot, aes(Dim1, Dim2)) + 
      geom_point(aes(colour = Group)) + 
      geom_text_repel(aes(label=Sample), hjust=0, vjust=0) +
      theme_bw3() + 
      theme(aspect.ratio=1) +
      xlab(paste0("Leading logFC dim 1 (", round(mds$var.explained[1] * 10000)/100, "%)")) + 
      ylab(paste0("Leading logFC dim 2 (", round(mds$var.explained[2] * 10000)/100, "%)")) +
      scale_color_manual(values=group_colors)
    ggsave(paste0(prefix, ".mds.png"), g, width=6, height=5, units="in", dpi=300, bg="white")

    #https://support.bioconductor.org/p/133907/#133920
    # To make a PCA plot, simply use
    # plotMDS(x, gene.selection="common")
    pca <- plotMDS(dge, gene.selection="common", plot=FALSE)
    toplot <- data.frame(Dim1 = pca$x, Dim2 = pca$y, Sample = design_data$Sample, Group = design_data$DisplayGroup)
    g<-ggplot(toplot, aes(Dim1, Dim2)) + 
      geom_point(aes(colour = Group)) + 
      geom_text_repel(aes(label=Sample), hjust=0, vjust=0) +
      theme_bw3() + 
      theme(aspect.ratio=1) +
      xlab(paste0("PC1 (", round(pca$var.explained[1] * 10000)/100, "%)")) + 
      ylab(paste0("PC2 (", round(pca$var.explained[2] * 10000)/100, "%)")) +
      scale_color_manual(values=group_colors)
    ggsave(paste0(prefix, ".pca.png"), g, width=6, height=5, units="in", dpi=300, bg="white")

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

  save_volcano_plot(edgeR_out_table=out$table, 
                                  prefix=prefix, 
                                  useRawPvalue=useRawPvalue, 
                                  pvalue=pvalue, 
                                  foldChange=foldChange, 
                                  comparisonTitle=paste0(cellType, " : ", comp))

  cat("  cpm", "\n")
  log_cpm <- cpm(dge, log=TRUE)
  cpm_file = paste0(prefix, ".log_cpm.csv")
  write.csv(log_cpm, file=cpm_file, quote=F)

  cat("  heatmap", "\n")
  topVarGenes <- head(order(rowVars(log_cpm), decreasing = TRUE), min(2000, nrow(log_cpm)))
  top_mat  <- log_cpm[ topVarGenes, ]

  suppressPackageStartupMessages(library("ComplexHeatmap"))

  drawHCA(paste0(prefix, ".log_cpm"), top_mat, design_data, group_colors, usePearsonInHCA=FALSE)
  
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

  gseaFile = write_gsea_rnk_by_loose_criteria(dge_all, groups, design, prefix)

  curDF<-data.frame("prefix"=prefix, "cellType"=cellType, "comparison"=comp, "betweenCluster"=0, "sampleInGroup"=0, "deFile"=dge_filename, "sigFile"=sigFile, "sigGenenameFile"=sigGenenameFile, "gseaFile"=gseaFile, "designFile"=design_file, "cpmFile"=cpm_file)
  if(is.null(result)){
    result<-curDF
  }else{
    result<-rbind(result, curDF)
  }
}

write.csv(result, file=paste0(outFile, ".edgeR.files.csv"), quote=F)



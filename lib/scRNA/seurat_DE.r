rm(list=ls()) 
outFile='P10940'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/nd_seurat_sct2_merge_dr0.2_3_choose/result/P10940.final.rds'
parFile2='/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/nd_seurat_sct2_merge_dr0.2_3_choose/result/P10940.meta.rds'
parFile3=''


setwd('/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/nd_seurat_sct2_merge_dr0.2_3_choose_SeuratDE_DESeq2_bySample_inCluster/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
source("scRNA_func.r")
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
DE_by_cell=is_one(myoptions$DE_by_cell)
DE_method=myoptions$DE_method

if(!DE_by_cell) {
  min_cell_per_sample=as.numeric(myoptions$filter_min_cell_per_sample)
}

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
  cat("Reading object from", parFile1, "\n")
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

detail_folder = paste0(outFile, ".details/")
if(!dir.exists(detail_folder)){
  dir.create(detail_folder)
}

drawHCA<-function(prefix, top_mat, designdata, group_colors, usePearsonInHCA=FALSE){
  mat_scaled = t(scale(t(top_mat)))

  ha=HeatmapAnnotation( Group=designdata$DisplayGroup,
                        col=list(Group=group_colors),
                        annotation_legend_param = list(Group = list(ncol = 1, title = "Group", title_position = "topleft")))

  if(usePearsonInHCA){
    clustering_distance_columns="pearson"
  }else{
    clustering_distance_columns="euclidean"
  }

  ignored = draw_heatmap_png( filepath=paste0(prefix, ".heatmap.png"), 
                              htdata=mat_scaled, 
                              name="zscore", 
                              show_row_names=FALSE, 
                              show_column_names=TRUE,
                              save_rds=FALSE,
                              clustering_distance_columns=clustering_distance_columns,
                              top_annotation=ha,
                              show_row_dend=FALSE)
}

detail_prefix = paste0(detail_folder, outFile)

meta<-obj@meta.data
mt<-data.frame(table(meta[,cluster_name], meta$orig.ident))
colnames(mt)<-c("cell_type", "sample", "num_cell")
write.csv(mt, paste0(detail_prefix, ".num_cell.csv"), row.names=F)

cts=sort(unique(meta[,cluster_name]))
cat("get design matrix ...\n")
designMatrix<-NULL

comp_dic=list()
comp <-comparisonNames[1]
for (comp in comparisonNames){
  comp_groups<-comparisons[comparisons$Comparison==comp,]
  comp_options = split(comp_groups$Value, comp_groups$Key)
  
  if("groups" %in% names(comp_options)){
    if(!is.null(myoptions$group_column) & (myoptions$group_column != "")){
      if(!myoptions$group_column %in% colnames(obj@meta.data)){
        stop(paste0("group_column ", myoptions$group_column, " not found in meta.data of ", parFile1))
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
  comp_dic[[comp]] = list(control_names=control_names, sample_names=sample_names, controlGroup=controlGroup, sampleGroup=sampleGroup, covariances=covariances)
}

result<-NULL

idx<-2
for (idx in c(1:length(cts))){
  ct = cts[idx]
  cts_file = celltype_to_filename(ct)

  #get the cells in the cell type
  ct_cells = rownames(meta)[meta[,cluster_name] == ct]
  ct_obj = subset(obj, cells = ct_cells)

  for (comp in comparisonNames){
    prefix = paste0(detail_folder, cts_file, ".", comp)
    cat("---", prefix, "---\n")

    cur_lst = comp_dic[[comp]]

    control_names = cur_lst$control_names
    sample_names = cur_lst$sample_names
    controlGroup = cur_lst$controlGroup
    sampleGroup = cur_lst$sampleGroup

    #get the cells in the comparison 
    comp_cells=colnames(ct_obj)[ct_obj$orig.ident %in% c(control_names, sample_names)]
    comp_obj=subset(ct_obj, cells=comp_cells)

    if(!DE_by_cell){
      old_samples<-unique(comp_obj$orig.ident)
      cat("  AggregateExpression for pseudo bulk DE analysis ...\n")
      comp_obj <- AggregateExpression(comp_obj, assays = "RNA", return.seurat = T, group.by = c("orig.ident"))
      if(!all(colnames(comp_obj) %in% old_samples)){
        if(all(gsub("-", "_", colnames(comp_obj)) %in% old_samples)){
          colnames(comp_obj)<-gsub("-", "_", colnames(comp_obj))
          comp_obj@meta.data$orig.ident<-gsub("-", "_", comp_obj@meta.data$orig.ident)
        }else{
          stop("AggregateExpression changed the sample names. Old samples: ", paste0(old_samples, collapse=","), " New samples: ", paste0(colnames(comp_obj), collapse=","))
        }
      }
    }

    comp_obj@meta.data$Group=factor(ifelse(comp_obj$orig.ident %in% control_names, "control", "sample"), levels=c("control", "sample"))
    comp_obj@meta.data$DisplayGroup=factor(ifelse(comp_obj$orig.ident %in% control_names, controlGroup, sampleGroup), levels=c(controlGroup, sampleGroup))

    if(!is.null(covariances)){
      for(cov_name in covariances){
        comp_obj@meta.data[,cov_name] = unlist(covariances_tbl[comp_obj@meta.data$orig.ident, cov_name])
      }
    }

    designdata=unique(comp_obj@meta.data[,c("orig.ident", "Group", "DisplayGroup", covariances)]) |>
      dplyr::rename(Sample=orig.ident)
    rownames(designdata)<-designdata$Sample
    designdata = designdata[c(control_names, sample_names),]

    design_file<-paste0(prefix, ".design")
    write.csv(designdata, file=design_file, row.names=F, quote=F)
    
    group_colors=c("blue", "red")
    names(group_colors)<-c(controlGroup, sampleGroup)

    sample_colors=group_colors[designdata$DisplayGroup]
    names(sample_colors)=designdata$Sample

    if(!DE_by_cell){
      cat(" draw figures using edgeR pesudo counts ...\n")
      library(edgeR)
      counts<-GetAssayData(comp_obj, assay="RNA", layer="counts")
      counts<-counts[,designdata$Sample]
      counts_file = paste0(prefix, ".counts.csv")
      write.csv(counts, file=counts_file, quote=F)

      groups<-designdata$DisplayGroup
      dge <- DGEList(counts=counts, group=groups)
    
      cat("  calcNormFactors", "\n")
      dge<-calcNormFactors(dge,method = "TMM")

      cat("  plotMDS", "\n")
      mds <- plotMDS(dge, plot=FALSE)
      toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Sample = designdata$Sample, Group = designdata$DisplayGroup)
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
      toplot <- data.frame(Dim1 = pca$x, Dim2 = pca$y, Sample = designdata$Sample, Group = designdata$DisplayGroup)
      g<-ggplot(toplot, aes(Dim1, Dim2)) + 
        geom_point(aes(colour = Group)) + 
        geom_text_repel(aes(label=Sample), hjust=0, vjust=0) +
        theme_bw3() + 
        theme(aspect.ratio=1) +
        xlab(paste0("PC1 (", round(pca$var.explained[1] * 10000)/100, "%)")) + 
        ylab(paste0("PC2 (", round(pca$var.explained[2] * 10000)/100, "%)")) +
        scale_color_manual(values=group_colors)
      ggsave(paste0(prefix, ".pca.png"), g, width=6, height=5, units="in", dpi=300, bg="white")

      cat("  cpm", "\n")
      log_cpm <- cpm(dge, log=TRUE)
      cpm_file = paste0(prefix, ".log_cpm.csv")
      write.csv(log_cpm, file=cpm_file, quote=F)

      cat("  heatmap", "\n")
      topVarGenes <- head(order(rowVars(log_cpm), decreasing = TRUE), min(2000, nrow(log_cpm)))
      top_mat  <- log_cpm[ topVarGenes, ]

      suppressPackageStartupMessages(library("ComplexHeatmap"))

      drawHCA(paste0(prefix, ".log_cpm"), top_mat, designdata, group_colors, usePearsonInHCA=FALSE)
    }else{
      log_cpm = NULL
      cpm_file = NULL
    }

    curdf<-data.frame(prefix=prefix, cellType=ct, comparison=comp, design=design_file, stringsAsFactors = F)
    if (is.null(designMatrix)){
      designMatrix = curdf
    }else{
      designMatrix = rbind(designMatrix, curdf)
    }

    Idents(comp_obj) <- "Group"

    out_table = FindMarkers(object = comp_obj, 
                      ident.1 = "sample", 
                      ident.2 = "control",
                      test.use = "DESeq2")
    out_table = out_table[!is.na(out_table$p_val),,drop=FALSE]

    #make it consistant with edgeR result
    out_table = out_table |>
                dplyr::rename(PValue=p_val, 
                              FDR=p_val_adj, 
                              logFC=avg_log2FC)


    dge_filename <-paste0(prefix, ".csv")
    write.csv(out_table, file=dge_filename, quote=F)

    if(useRawPvalue){
      sigout<-out_table[(out_table$PValue<=pvalue) & (abs(out_table$logFC)>=log2(foldChange)),]
    }else{
      sigout<-out_table[(out_table$FDR<=pvalue) & (abs(out_table$logFC)>=log2(foldChange)),]
    }
    sigFile<-paste0(prefix, ".sig.csv")
    write.csv(sigout, file=sigFile, quote=F)
    
    sig_genes<-data.frame(gene=rownames(sigout), stringsAsFactors = F)
    sigGenenameFile<-paste0(prefix, ".sig_genename.txt")
    write.table(sig_genes, file=sigGenenameFile, row.names=F, col.names=F, sep="\t", quote=F)

    rankout<-data.frame(gene=rownames(out_table), sigfvalue=sign(out_table$logFC) * (-log10(out_table$PValue)))
    rankout<-rankout[order(rankout$sigfvalue, decreasing=TRUE),]
    gseaFile<-paste0(prefix, "_GSEA.rnk")
    write.table(rankout, file=gseaFile, row.names=F, col.names=F, sep="\t", quote=F)

    save_volcano_plot(edgeR_out_table=out_table, 
                      prefix=prefix, 
                      useRawPvalue=useRawPvalue, 
                      pvalue=pvalue, 
                      foldChange=foldChange, 
                      comparisonTitle=paste0(ct, " : ", comp))

    if(nrow(sig_genes) > 0){
      sig_gene = sig_genes$gene[1]

      cell_obj=subset(ct_obj, cells=comp_cells)

      g<-get_sig_gene_figure(cell_obj, out_table, designdata, sig_gene, DE_by_cell=DE_by_cell, is_between_cluster=FALSE, log_cpm=log_cpm)
      png(paste0(prefix, ".top_1_gene.png"), width=3000, height=2500, res=300)
      print(g)
      dev.off()
    }

    curDF<-data.frame("prefix"=prefix, "cellType"=ct, "comparison"=comp, "betweenCluster"=0, "sampleInGroup"=0, "deFile"=dge_filename, "sigFile"=sigFile, "sigGenenameFile"=sigGenenameFile, "gseaFile"=gseaFile, "designFile"=design_file, "cpmFile"=cpm_file)
    if(is.null(result)){
      result<-curDF
    }else{
      result<-rbind(result, curDF)
    }
  }
}

write.csv(result, file=paste0(outFile, ".SeuratDE.files.csv"), quote=F)

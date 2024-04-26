rm(list=ls()) 
outFile='SADIE_adipose'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/nobackup/shah_lab/shengq2/20240304_mona_scRNA_SADIE/20240422_composition/1_prepare_data/adipose_myeloid.rds'
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20240304_mona_scRNA_SADIE/20240422_composition/2_dcats/result')

### Parameter setting end ###

source("scRNA_func.r")
#########################
#### 2024.03.06 Jing Yang
#### using DCATS with similarity matrix
##########################

library(Seurat)
#BiocManager::install('holab-hku/DCATS')
library(DCATS)

opt=fread("fileList1.txt",header=F)
myoptions=split(opt$V1, opt$V2)
myoptions$by_sctransform=is_one(myoptions$by_sctransform)
myassay=ifelse(myoptions$by_sctransform, "SCT", "RNA")
reduction=myoptions$DCATS_reduction
sample_column=myoptions$DCATS_sample_column
celltype_column=myoptions$DCATS_celltype_column
reference_celltype=myoptions$DCATS_reference_celltype

obj <- readRDS(parFile1)
ndim=ncol(Embeddings(obj, reduction = reduction))

obj <- FindNeighbors(obj,reduction=reduction,dims=1:ndim,assay=myassay)

composition_data<-FetchData(obj, vars=c(sample_column, celltype_column))
composition<-as.matrix(table(composition_data[,sample_column], composition_data[,celltype_column]))
write.csv(composition, file=paste0(outFile,".ct_sample.csv"), row.names = T)

if(reference_celltype != ""){
  if(!reference_celltype %in% colnames(composition)){
    stop(paste0(reference_celltype, " must be in object column ", celltype_column))
  }
}

groups<-fread("fileList2.txt",header=F, data.table=F) %>% 
  dplyr::rename("Sample"=V1, "Group"=V2)

comparison_data<-fread("fileList3.txt", header=F, data.table=FALSE) %>% 
  dplyr::filter(V2=="groups") %>%
  dplyr::select(V1, V3) %>%
  dplyr::rename("Group"=V1, "Comparison"=V3)

comparisons <- unique(comparison_data$Comparison)

results <- list()
k=1
for (k in 1:length(comparisons)) {
  cur_comp = comparisons[k]
  print(cur_comp)

  comp_groups=comparison_data %>% dplyr::filter(Comparison==cur_comp)
  comp_samples=groups[groups$Group %in% comp_groups$Group,]
  rownames(comp_samples)=comp_samples$Sample
  write.csv(comp_samples, file=paste0(outFile, ".", cur_comp, ".design.csv"), row.names = F)

  count_mat <- composition[comp_samples$Sample,]
  write.csv(t(count_mat), file=paste0(outFile, ".", cur_comp, ".cell_count.csv"), row.names = T)

  perc_mat=t(count_mat)
  perc_mat=round(t(t(perc_mat)/colSums(perc_mat)) * 100, 2)
  write.csv(perc_mat, file=paste0(outFile, ".", cur_comp, ".cell_perc.csv"), row.names = T)

  ## Percentage Boxplot
  box_width=7
  box_height=5
  mtbl=reshape2::melt(as.matrix(perc_mat))
  mtbl$Condition=comp_samples[mtbl$Var2,"Group"]
  mtbl=mtbl[order(mtbl$value, decreasing=TRUE),]
  mtbl$Var1=factor(mtbl$Var1, levels=rownames(perc_mat))
  g<-ggplot(mtbl, aes(x=Var1,y=value,color=Condition)) + 
    geom_boxplot() + 
    theme_bw3() + 
    xlab(NULL) + 
    ylab("Percentage of cell") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  cell_png=paste0(outFile, ".", cur_comp, ".cell_perc.png")
  ggsave(cell_png, g, width=box_width, height=box_height, dpi=300, units="in", bg="white")

  design_mat <- data.frame(condition = comp_samples$Group)
  cells =rownames(obj@meta.data)[obj@meta.data[,sample_column] %in% rownames(count_mat)]
  obj_sub <- subset(obj, cells=cells)
  
  print('knn_simMat')
  if(myoptions$by_sctransform){
    knn_mat <- knn_simMat(obj_sub@graphs$SCT_snn,obj_sub@meta.data[,celltype_column])
  }else{
    knn_mat <- knn_simMat(obj_sub@graphs$RNA_snn,obj_sub@meta.data[,celltype_column])
  }
  
  print('dcats_GLM')
  cur_result = dcats_GLM( count_mat=count_mat, 
                          design_mat=design_mat, 
                          similarity_mat=knn_mat, 
                          reference=reference_celltype)
  cur_tbl=as.data.frame(cur_result)
  colnames(cur_tbl)=names(cur_result)
  cur_tbl$comparison=cur_comp
  cur_tbl$cluster=rownames(cur_tbl)
  cur_tbl=cur_tbl %>% dplyr::select(cluster, comparison, everything())

  write.csv(cur_tbl, file=paste0(outFile,".",cur_comp,".DCATS.csv"), row.names = F)

  results[[comparisons[k]]] <- cur_tbl
}
saveRDS(results,file=paste0(outFile,".DCATS.rds"))
combined=do.call(rbind, results)
write.csv(combined, file=paste0(outFile,".DCATS.csv"), row.names = F)


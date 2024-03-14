rm(list=ls()) 
outFile='int_papaer_crs'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/nobackup/h_turner_lab/yangj22/20231031_integrate_a_paper_and_20230427_7114_8822_scRNA_hg38_vst2/result20231214/seurat_sct2_merge_dr0.5_3_choose/result/int_papaer_crs.final.rds'
parFile2=''
parFile3=''


setwd('/nobackup/h_turner_lab/yangj22/20231031_integrate_a_paper_and_20230427_7114_8822_scRNA_hg38_vst2/result20231214/seurat_sct2_merge_dr0.5_3_choose_DCATS/result')

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

obj <- readRDS(parFile1)
obj <- FindNeighbors(obj,reduction="pca",dims=1:30,assay=myassay)

composition_data<-FetchData(obj, vars=c("orig.ident", myoptions$celltype_column))
composition<-as.matrix(table(composition_data$orig.ident, composition_data[,myoptions$celltype_column]))
write.csv(composition, file=paste0(outFile,".ct_sample.csv"), row.names = T)

reference_celltype=NULL
if(myoptions$reference_celltype != ""){
  if(!myoptions$reference_celltype %in% colnames(composition)){
    stop(paste0(myoptions$reference_celltype, " must be in object column ", myoptions$celltype_column))
  }
  reference_celltype=myoptions$reference_celltype
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
  mtbl=melt(as.matrix(perc_mat))
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
  obj_sub <- subset(obj, orig.ident %in% rownames(count_mat))
  
  print('knn_simMat')
  if(myoptions$by_sctransform){
    knn_mat <- knn_simMat(obj_sub@graphs$SCT_snn,obj_sub@meta.data[,myoptions$celltype_column])
  }else{
    knn_mat <- knn_simMat(obj_sub@graphs$RNA_snn,obj_sub@meta.data[,myoptions$celltype_column])
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

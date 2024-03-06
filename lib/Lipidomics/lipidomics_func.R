library(data.table)
library(dplyr)

make_valid_name<-function(old_name){
  new_name <- gsub(" ", "_", old_name)
  new_name <- gsub("-", "_", new_name)
  new_name <- gsub("\\+", "_", new_name)
  is_number <- grepl("^[0-9]", new_name)
  new_name[is_number] <- paste0("X", new_name[is_number])
  return(new_name)
}

read_lipid_tsv_and_filter<-function(filename, filePrefix, qc_group="QC", blank_groups=c("Solvent", "Blank"), not_sample_groups=c("AIS", "dBlank")){
  modePrefix=".*[POS|NEG]_"

  header<-data.frame(t(fread(filename, nrows=5, data.table=FALSE)))

  nameCol=min(which("Metabolite name"==header$X5))
  ontologyCol=which("Ontology"==header$X5)
  classCol=which("Class"==header$X1)
  aveCols=which("Average"==header$X4)
  stdevCols=which("Stdev"==header$X4)
  sampleCols=(classCol + 1):(min(aveCols) - 1)
  sampleStartCol=min(sampleCols)

  colnames(header)[1:4]<-header[sampleStartCol-1,c(1:4)]
  
  meta=header[!is.na(header$Class) & !(header$Class %in% c("", "Class")),]
  colnames(meta)[5]="Sample"
  meta$Sample=make.names(gsub(modePrefix, "", meta$Sample))
  meta$Class=make.names(meta$Class)
  meta <- meta %>% dplyr::select("Sample", everything())

  write.csv(meta, paste0(filePrefix, ".meta.csv"), row.names=FALSE)

  all_not_sample_groups=c(qc_group, blank_groups, not_sample_groups)
  sample_meta=meta[!(meta[, "Class"] %in% all_not_sample_groups) & meta[, "File type"] == "Sample",]

  sample_meta$Class_Sample <- make_valid_name(paste0(sample_meta$Class, "_", sample_meta$Sample))
  
  write.csv(sample_meta, paste0(filePrefix, ".sample_meta.csv"), row.names=FALSE)

  df<-fread(filename, skip=4, data.table=FALSE)
  #convert all integer64 to numeric, otherwise, calculation of mean and sd will be wrong, although we will not calculate it here
  for(i in 1:ncol(df)){
    if(class(df[,i])=="integer64"){
      df[,i]=as.numeric(df[,i])
    }
  }

  colnames(df)[nameCol]<-"Metabolite.name"

  #remove "POS_" or "NEG_" prefix from sample names to make sure the final sample names in both POS and NEG file would be idetical
  colnames(df)[sampleCols]<-gsub(modePrefix,"",colnames(df)[sampleCols])

  res=df[,c(nameCol, ontologyCol, which(colnames(df) %in% sample_meta$Sample))]

  qc_samples=meta[meta[,"Class"]==qc_group, "Sample"]
  if(length(qc_samples)==0){
    stop(paste0("No QC sample found for ", qc_group))
  }
  if(length(qc_samples)==1){
    msg=paste0("Only one QC sample found for ", qc_group, " group: ", qc_samples)
    cat("\n<mark>", msg, "</mark>\n")
  }

  blank_samples=meta[meta[,"Class"] %in% blank_groups, "Sample"]
  if(length(blank_samples)==0){
    stop(paste0("No Blank sample found for ", blank_groups))
  }
  if(length(blank_samples)==1){
    msg=paste0("Only one blank sample found for ", paste0(blank_groups, collapse = "/"), " group: ", blank_samples)
    cat("\n<mark>", msg, "</mark>\n")
  }

  if(FALSE){
    #we would not calculate the QC and Blank here, because we would like to use the original data
    res$mean_QC=apply(df[,qc_samples], 1, mean)
    res$mean_Blank=apply(df[,blank_samples], 1, mean)
    res$stdev_QC=apply(df[,qc_samples], 1, sd)
  }else{
    aver_df=df[,aveCols]
    stdev_df=df[,stdevCols]

    res$mean_QC=aver_df[,qc_group]
    
    blank_index=which(colnames(aver_df) %in% blank_groups)
    res$mean_Blank=aver_df[,blank_index]
    
    res$stdev_QC=stdev_df[,qc_group]
  }
  res$RSD_QC=100 * res$stdev_QC / res$mean_QC
  res$blank_QC_ratio=res$mean_Blank / res$mean_QC

  write.csv(res, paste0(filePrefix, ".with_QC.csv"), row.names=FALSE)

  filter_tb=data.frame(matrix(ncol=2, nrow=8))
  colnames(filter_tb)=c("Filter", "Count")
  filter_tb[1,]=c("Total features", nrow(res))  

  res <- res %>% filter(Metabolite.name != 'Unknown')
  filter_tb[2,]=c('remove "Unknown"', nrow(res))  

  res <- res %>% filter(!grepl('w/o', Metabolite.name))
  filter_tb[3,]=c('remove "w/o MS2"', nrow(res))  

  res <- res %>% filter(Ontology != 'Others') %>% filter(!grepl('RIKEN', Metabolite.name))
  filter_tb[4,]=c('remove other unidentified features', nrow(res))  

  res <- res %>%  filter(!grepl('\\(d7\\)', Metabolite.name)) %>%
                  filter(!grepl('\\-d7\\s', Metabolite.name)) %>%
                  filter(!grepl('\\(d9\\)', Metabolite.name)) %>%
                  filter(!grepl('\\-d9\\s', Metabolite.name)) 

  filter_tb[5,]=c('remove d7/d9', nrow(res)) 

  res=res[res$mean_QC > 0,]
  filter_tb[6,]=c("mean(QC) > 0", nrow(res))  

  res=res[res$RSD_QC < 25,]
  filter_tb[7,]=c("RSD(QC) < 25%", nrow(res))  

  res=res[res$blank_QC_ratio < 0.1,]
  filter_tb[8,]=c("mean(Blank)/mean(QC) < 0.1", nrow(res))  
  write.csv(filter_tb, paste0(filePrefix, ".filter.csv"), row.names=FALSE)

  dupnames=res$Metabolite.name[duplicated(res$Metabolite.name)]
  dup_tb=res[res$Metabolite.name %in% dupnames,]
  write.csv(dup_tb, paste0(filePrefix, ".dup.csv"), row.names=FALSE)

  dup_df=df[df$Metabolite.name %in% dupnames,]
  write.csv(dup_df, paste0(filePrefix, ".dup_all.csv"), row.names=FALSE)

  res = res %>% dplyr::select(-c(mean_QC, mean_Blank, stdev_QC, RSD_QC, blank_QC_ratio))

  return(list(res=res, all_meta=meta, sample_meta=sample_meta, filter_tb=filter_tb))
}

draw_pheatmap <- function(data_df, groups, annotation_colors, title, file_prefix, heatmap_width, heatmap_height, annotation_legend=FALSE, show_rownames=nrow(data_df) <= 30){
  cellheight = ifelse(show_rownames, 12, NA)
  max_height=ifelse(show_rownames, 7, 6)
  heatmap_height = ifelse(is.na(heatmap_height), max(4, min(max_height, nrow(data_df) * 0.2 + 3)), heatmap_height)
  png_file = paste0(file_prefix, ".heatmap.png")
  pheatmap( data_df, 
            main = title,
            scale = "row", 
            show_rownames = show_rownames, 
            annotation_col = groups, 
            annotation_colors = annotation_colors, 
            annotation_legend = annotation_legend, 
            annotation_names_col = F, 
            border_color = NA,  
            fontsize = 10, 
            cellheight = cellheight, 
            cluster_cols = T,
            filename=png_file,
            width=heatmap_width,
            height=heatmap_height)
  return(png_file)
}

draw_pca <- function(data, groups, title, file_prefix, pca_width, pca_height, colors=NULL){
  data_t <- t(data)
  pca <- prcomp(data_t)
  pcadata<-data.frame(pca$x)
  pcadata=as.data.frame(scale(pcadata))
  sum_pca <- summary(pca)$importance
  pcadata$groups <- groups
  
  p <- ggplot(pcadata, aes(x=PC1, y=PC2, color=groups)) + #, label = rownames(pcadata))) 
    geom_point(size = 4) + 
    theme_bw() +
    labs(x=paste0("PCA1 (", round(sum_pca[2,1] * 100, digits = 1), "%)"),
         y=paste0("PCA2 (", round(sum_pca[2,2] * 100, digits = 1), "%)"),
         title=title) +
    theme(axis.text = element_text(size = 16), 
          axis.title = element_text(size=18, face = "bold"), 
          legend.text=element_text(size=12), 
          legend.title= element_blank(),
          plot.title = element_text(hjust = 0.5),
          aspect.ratio=1)
  
  if(length(unique(groups)) < 10 & !is.null(colors)){
    p <- p + #scale_color_brewer(palette = "Set1")
    scale_colour_manual(values = c(colors))
  }

  png_file = paste0(file_prefix, ".pca.png")
  ggsave(png_file, plot = p, width = pca_width, height = pca_height, dpi = 300, units = "in", bg="white")
  return(png_file)
}

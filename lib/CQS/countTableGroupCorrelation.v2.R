rm(list=ls()) 
outFile=''
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1=''
parFile2=''
parFile3=''
outputPdf<-FALSE;outputPng<-TRUE;outputTIFF<-FALSE;showVolcanoLegend<-FALSE;usePearsonInHCA<-TRUE;showLabelInPCA<-FALSE;useGreenRedColorInHCA<-FALSE;top25cvInHCA<-FALSE;

setwd('/nobackup/vickers_lab/projects/20240206_Linton_11055_bulkRNA_human/genetable/result')

### Parameter setting end ###

source("countTableVisFunctions.R")

library(data.table)

options(bitmapType='cairo')

library(DESeq2)  
library(RColorBrewer)
library(colorRamps)
library(genefilter)
library(limma)
library(dplyr)
library(ggplot2)
library(tibble)
library(cowplot)
library(ComplexHeatmap)

if(exists("parSampleFile4")){
  myoptions_tbl<-read.table(parSampleFile4, sep="\t", stringsAsFactors = F)
  myoptions<-split(myoptions_tbl$V1, myoptions_tbl$V2)

  draw_all_groups_in_HCA<-is_one(myoptions$draw_all_groups_in_HCA)
  draw_umap<-is_one(myoptions$draw_umap)
  heatmap_cexCol = myoptions$heatmap_cexCol
  if(!is.na(heatmap_cexCol)){
    if(heatmap_cexCol == ""){
      heatmap_cexCol<-NA
    }else{
      heatmap_cexCol<-as.numeric(heatmap_cexCol)
    }
  }

  if("n_first" %in% names(myoptions)){
    n_first = as.numeric(myoptions$n_first)
  }else{
    n_first = -1
  }
}else{
  draw_all_groups_in_HCA<-FALSE
  draw_umap<-FALSE
  heatmap_cexCol<-NA
  n_first = -1
}

countTableFileList<-parSampleFile1
groupFileList<-parSampleFile2
colorFileList<-parSampleFile3

fixColorRange<-TRUE

geneFile<-parFile1
totalCountFile<-parFile3

covarianceFile<-ifelse(exists("parFile4"), parFile4, "")
if(file.exists((covarianceFile))){
  covariances<-fread(covarianceFile, header=T, data.table=F)
  has_batch<-"batch" %in% colnames(covariances)
  if(has_batch){
    batch_map<-unlist(split(covariances$batch, covariances$Sample))
  }
}else{
  has_batch<-FALSE
}

if(!exists("onlySamplesInGroup")){
  onlySamplesInGroup=TRUE
}

if(!exists("outputPdf")){
  outputPdf<-FALSE
}

if(!exists("outputPng") | !outputPdf ){
  outputPng<-TRUE
}

outputFormat<-c()
if(outputPdf){
  outputFormat<-c("PDF")
}
if(outputPng){
  outputFormat<-c(outputFormat, "PNG")
}

if(!exists("hasRowNames")){
  hasRowNames<-NA
}

if(exists("useGreenRedColorInHCA") && useGreenRedColorInHCA){
  hmcols <- colorRampPalette(c("green", "black", "red"))(256)
}else{
  hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(256)
}

if(!exists("useGroupAsBatch")){
  useGroupAsBatch<-FALSE
}

if(!exists("showLabelInPCA")){
  showLabelInPCA<-TRUE
}

if(!exists("transformTable")){
  transformTable<-FALSE
}

if(!exists("suffix")){
  suffix<-""
}

if(!exists("idIndex")){
  idIndex<-1
}

task_suffix<-suffix

#outputDirectory<-"."
output_include_folder_name<-1
if(!exists("outputDirectory")){
  outputDirectory<-""
}

if (!exists("minMedian")){
  minMedian<-1
}

if (!exists("minMedianInGroup")){
  minMedianInGroup<-1
}

#source("/home/zhaos/source/ngsperl/lib/CQS/countTableVisFunctions.R")

##Solving node stack overflow problem start###
#when there are too many genes, drawing dendrogram may failed due to node stack overflow,
#It could be solved by forcing stats:::plotNode to be run as interpreted code rather then byte-compiled code via a nasty hack.
#http://stackoverflow.com/questions/16559250/error-in-heatmap-2-gplots/25877485#25877485

# Convert a byte-compiled function to an interpreted-code function 
unByteCode <- function(fun)
{
  FUN <- eval(parse(text=deparse(fun)))
  environment(FUN) <- environment(fun)
  FUN
}

# Replace function definition inside of a locked environment **HACK** 
assignEdgewise <- function(name, env, value)
{
  unlockBinding(name, env=env)
  assign( name, envir=env, value=value)
  lockBinding(name, env=env)
  invisible(value)
}

# Replace byte-compiled function in a locked environment with an interpreted-code
# function
unByteCodeAssign <- function(fun)
{
  name <- gsub('^.*::+','', deparse(substitute(fun)))
  FUN <- unByteCode(fun)
  retval <- assignEdgewise(name=name,
                           env=environment(FUN),
                           value=FUN
  )
  invisible(retval)
}

# Use the above functions to convert stats:::plotNode to interpreted-code:
unByteCodeAssign(stats:::plotNode)

# Now raise the interpreted code recursion limit (you may need to adjust this,
#  decreasing if it uses to much memory, increasing if you get a recursion depth error ).
options(expressions=5e4)
drawPCA<-function(filename, rldmatrix, showLabelInPCA, groups, groupColors, outputFormat, width_inch=7, height_inch=6,scalePCs=TRUE){
  genecount<-nrow(rldmatrix)
  if(genecount > 2){
    pca<-prcomp(t(rldmatrix))
    supca<-summary(pca)$importance

    npc = min(10, ncol(supca))
    pca_res = t(supca)[c(1:npc),] %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column("PC") %>% 
      dplyr::rename("Proportion" = "Proportion of Variance", "Cumulative" = "Cumulative Proportion")

    pca_res$PC <- factor(pca_res$PC, levels = pca_res$PC)

    g1<-ggplot(pca_res, aes(x=PC, y=Proportion)) +
      geom_bar(stat="identity", fill="steelblue") +
      geom_text(aes(label=Proportion), vjust=-0.3, size=3.5) +
      scale_y_continuous(labels = scales::percent) +
      labs(x = "Principal Components", y = "Proportion of Variance", title = "Proportion of Variance Explained by Each PC") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    g2<-ggplot(pca_res, aes(x=PC, y=Cumulative)) + 
      geom_bar(stat="identity", fill="steelblue") +
      geom_text(aes(label=Cumulative), vjust=-0.3, size=3.5) +
      scale_y_continuous(labels = scales::percent) +
      labs(x = "Principal Components", y = "Cumulative Proportion of Variance", title = "Cumulative Proportion of Variance Explained by PCs") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

    g<-ggdraw() +
      draw_plot(g1, x = 0, y = .5, width = 1, height = .5) +
      draw_plot(g2, x = 0, y = 0, width = 1, height = .5)

    ggsave(paste0(filename, ".bar.png"), g, width = 8, height = 8, dpi = 300, units = "in", bg="white")

    pcadata<-data.frame(pca$x)
    if (scalePCs) {
      pcadata=as.data.frame(scale(pcadata))
    }
    pcalabs=paste0(colnames(pcadata), "(", round(supca[2,] * 100), "%)")
    pcadata["sample"]<-row.names(pcadata)
    if(!is.na(groups[1])){
      pcadata$group<-groups
    }
    
    if(showLabelInPCA){
      g <- ggplot(pcadata, aes(x=PC1, y=PC2, label=sample)) + 
        geom_text(vjust=-0.6, size=4)
    }else{
      g <- ggplot(pcadata, aes(x=PC1, y=PC2)) +
        theme(legend.position="top")
    }
    if(!is.na(groups[1])){
      g<-g+geom_point(aes(col=group), size=4) + 
        scale_colour_manual(name="",values = groupColors)
    }else{
      g<-g+geom_point(size=4) 
    }
    g<-g+scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) +
      scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) + 
      geom_hline(aes(yintercept=0), linewidth=.2) + 
      geom_vline(aes(xintercept=0), linewidth=.2) + 
      xlab(pcalabs[1]) + ylab(pcalabs[2]) +
      theme_bw2() + theme(aspect.ratio=1)

    save_ggplot2_plot(file_prefix=filename, 
      outputFormat=outputFormat, 
      width_inch=width_inch, 
      height_inch=height_inch, 
      g)
  }
}

##Solving node stack overflow problem end###

if (geneFile!="") { #visualization based on genes in geneFile only
  genes<-read.table(geneFile, sep="\t", header=F, stringsAsFactors=F)$V1
  cat("There are", length(genes), "genes in gene file.\n")
}else{
  genes<-NA
}

if(colorFileList != ""){
  colorsTab<-read.table(colorFileList, sep="\t", header=F, stringsAsFactors=F)
  groupColors<-colorsTab$V1
  names(groupColors)<-colorsTab$V2
}else{
  groupColors<-NA
}

#start work:
countTableFileAll<-read.delim(countTableFileList,header=F,as.is=T,check.names=F)

missed_count_tables = c()
missed_count_tables_file<-"missed_count_tables.txt"

succeed_file<-"correlation.succeed"
if(file.exists(succeed_file)){
  unlink(succeed_file)
}

i<-2
for (i in 1:nrow(countTableFileAll)) {
  countTableFile<-countTableFileAll[i,1]
  #countTableFile<-paste0("C:/projects", countTableFile)

  if(!file.exists(countTableFile)){
    missed_count_tables<-c(missed_count_tables, countTableFile)
    next
  }
  
  if(outputDirectory==""){
    outputFilePrefix=countTableFile
  }else{
    bname=basename(countTableFile)
    if(output_include_folder_name){
      dpath=dirname(countTableFile)
      dname=basename(dpath)
      while(dname=="result"){
        dpath=dirname(dpath)
        dname=basename(dpath)
      }
      outputFilePrefix=paste0(outputDirectory, "/", dname, "." , bname)
    }else{
      outputFilePrefix=paste0(outputDirectory, "/", bname)
    }
  }
  
  countTableTitle<-countTableFileAll[i,2]
  
  cat("Reading",countTableFile, "\n")
  
  if(n_first != -1){
    count<-data.frame(fread(countTableFile, nrows=n_first), row.names=idIndex,check.names=FALSE)
  }else{
    if (grepl(".csv$",countTableFile)) {
      count<-read.csv(countTableFile,header=T,row.names=idIndex,as.is=T,check.names=FALSE)
    } else {
      count<-read.delim(countTableFile,header=T,row.names=idIndex,as.is=T,check.names=FALSE)
    }
  }
  
  if(transformTable){
    count<-t(count)
  }
  
  if (nrow(count)==0) {
    next;
  }
  
  if (ncol(count)<2) {
    next;
  }
  
  if(!is.na(genes)){
    if("Feature_gene_name" %in% colnames(count)){
      curgenes<-c(genes, rownames(count)[count$Feature_gene_name %in% genes])
    }else{
      curgenes<-genes
    }
  }else{
    curgenes<-NA
  }
  
  count[is.na(count)]<-0
  
  colClass<-sapply(count, class)
  countNotNumIndex<-which(colnames(count) %in% c("Feature_length","Location","Sequence") | (colClass!="numeric" & colClass!="integer"))
  if (length(countNotNumIndex)==0) {
    countNotNumIndex<-0;
  } else {
    countNotNumIndex<-max(countNotNumIndex)
  } #Please note here we only use columns after the rightest Non-Number column as count data
  countNum<-count[,c((countNotNumIndex+1):ncol(count))]
  countNum<-round(countNum,0)
  
  #remove genes with total reads 0
  countNum<-countNum[which(rowSums(countNum,na.rm=T)>0),]
  #remove samples with total reads 0
  countNum<-countNum[,which(colSums(countNum,na.rm=T)>0)]
  
  if (groupFileList!="") {
    sampleToGroup<-getSampleInGroup(groupFileList, colnames(countNum), useLeastGroups,onlySamplesInGroup=onlySamplesInGroup)
    if (onlySamplesInGroup) {
      if(!any(colnames(countNum) %in% as.character(sampleToGroup$V1))){
        stop(paste0("No file in count table was found in sample/group file:", groupFileList))
      }
      countNum<-countNum[,which(colnames(countNum) %in% as.character(sampleToGroup$V1))]
      #remove genes with total reads 0
      countNum<-countNum[which(rowSums(countNum,na.rm=T)>0),]
    }
  }else{
    sampleToGroup<-data.frame("V1"=colnames(countNum), "V2"="all", "V3"="all")
  }

  titles<-unique(sampleToGroup$V3)
  if(draw_all_groups_in_HCA){
    allSamples<-sort(unique(sampleToGroup$V1))
    conditionColors<-data.frame(Sample=allSamples)
    rownames(conditionColors)=allSamples

    title<-titles[1]
    for(title in titles){
      validSampleToGroup<-sampleToGroup[sampleToGroup$V3 == title,,drop=F]
      validSampleToGroup$V2<-factor(validSampleToGroup$V2, levels=sort(unique(validSampleToGroup$V2)))
      rownames(validSampleToGroup)<-validSampleToGroup$V1
      validSampleToGroup<-validSampleToGroup[allSamples,]
      
      if(title == "all"){
        title="Group"
      }
      
      groups<-validSampleToGroup$V2
      if (is.na(groupColors)){
        colors<-makeColors(length(unique(groups)))
      }else{
        colors<-groupColors
      }
      curColors<-data.frame(Group=colors[groups])
      conditionColors[,title]<-curColors
    }
    conditionColors<-conditionColors[,c(2:ncol(conditionColors)),drop=F]
    conditionColors<-as.matrix(conditionColors)
  }

  title<-titles[1]
  for(title in titles){
    cat("Processing", title, "samples\n")
    validSampleToGroup<-sampleToGroup[sampleToGroup$V3 == title,,drop=F]
    validSampleToGroup$V2<-factor(validSampleToGroup$V2, levels=sort(unique(validSampleToGroup$V2)))
    
    if(title != "all"){
      curSuffix = paste0(suffix, ".", title)
    }else{
      curSuffix = suffix
    }
    curSuffix<-gsub("\\s+", "_", curSuffix)
    
    write.table(validSampleToGroup, paste0(outputFilePrefix,curSuffix,".correlation.groups"),col.names=F, row.names=F, quote=F, sep="\t")
    
    #filter reads/genes by parameter
    validCountNum<-filterCountTable(countNum,validSampleToGroup,minMedian=minMedian,minMedianInGroup=minMedianInGroup)
    
    #normlize by total reads or VSD
    bNormalizeByCount=FALSE
    if(totalCountFile != ""){
      bNormalizeByCount=totalCountKey != "None"
    }

    if (bNormalizeByCount) { #normlize with total count *10^6
      totalCount<-read.csv(totalCountFile,header=T,as.is=T,row.names=1,check.names=FALSE)
      if(!(totalCountKey %in% rownames(totalCount))){
        if(!file.exists(parFile2)){
          stop(paste0(totalCountKey, " not exists in file ", totalCountFile))
        }
        totalCount<-read.csv(parFile2,header=T,as.is=T,row.names=1,check.names=FALSE)
        if(!(totalCountKey %in% rownames(totalCount))){
          stop(paste0(totalCountKey, " not exists in file ", totalCountFile, " and ", parFile2))
        }
      }
      totalCount<-unlist(totalCount[totalCountKey,])
      notValidSamples = colnames(validCountNum)[!(colnames(validCountNum) %in% names(totalCount))]
      if (length(notValidSamples) > 0){
        stop(paste0("Sample ", notValidSamples, " not found in total count table ", totalCountFile, "\n"))
      }
      countNumVsd<-10^6*t(t(validCountNum)/totalCount[colnames(validCountNum)])
      write.table(countNumVsd, paste0(outputFilePrefix,curSuffix,".RPM.txt"),col.names=NA, quote=F, sep="\t")
      countNumVsd<-log2(countNumVsd+1)
      ylab<-"log2(Mapped Reads per Million)"
    } else {
      dds=DESeqDataSetFromMatrix(countData = validCountNum, colData = as.data.frame(rep(1,ncol(validCountNum))),design = ~1)
      dds<-try(myEstimateSizeFactors(dds))
      vsdres<-try(dds<-DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
      if (class(vsdres) == "try-error") {
        message=paste0("Warning: varianceStabilizingTransformation function failed.\n",as.character(vsdres))
        warning(message)
        writeLines(message,paste0(outputFilePrefix,curSuffix,".vsd.warning"))
        next;
      }
      
      if(has_batch){
        saveRDS(dds, paste0(outputFilePrefix,curSuffix,".before_removeBatchEffect.vsd.rds"))

        countNumVsd<-assay(dds)
        colnames(countNumVsd)<-colnames(validCountNum)
        write.table(countNumVsd, paste0(outputFilePrefix,curSuffix,".before_removeBatchEffect.vsd.txt"),col.names=NA, quote=F, sep="\t")

        dds$batch<-batch_map[colnames(dds)]
        ntop = 3000

        png(paste0(outputFilePrefix,curSuffix,".before_removeBatchEffect.batch.png"), width=2000, height=2000, res=300)
        g<-plotPCA(dds, "batch", ntop=ntop) + theme_bw3()
        print(g)
        dev.off()

        hasMultipleGroup<-length(unique(validSampleToGroup$V2)) > 1
        if(hasMultipleGroup){
          stopifnot(all(rownames(colData(dds)) == validSampleToGroup$V1))
          dds$group = validSampleToGroup$V2

          png(paste0(outputFilePrefix,curSuffix,".before_removeBatchEffect.group.png"), width=2000, height=2000, res=300)
          g<-plotPCA(dds, "group", ntop=ntop) + theme_bw3()
          print(g)
          dev.off()

          png(paste0(outputFilePrefix,curSuffix,".before_removeBatchEffect.batch_group.png"), width=2000, height=2000, res=300)
          g<-plotPCA(dds, c("batch", "group"), ntop=ntop) + theme_bw3()
          print(g)
          dev.off()
        }

        assay(dds) <- limma::removeBatchEffect(assay(dds), dds$batch)
        png(paste0(outputFilePrefix,curSuffix,".after_removeBatchEffect.batch.png"), width=2000, height=2000, res=300)
        g<-plotPCA(dds, "batch", ntop=ntop) + theme_bw3()
        print(g)
        dev.off()

        if(hasMultipleGroup){
          png(paste0(outputFilePrefix,curSuffix,".after_removeBatchEffect.group.png"), width=2000, height=2000, res=300)
          g<-plotPCA(dds, "group", ntop=ntop) + theme_bw3()
          print(g)
          dev.off()

          png(paste0(outputFilePrefix,curSuffix,".after_removeBatchEffect.batch_group.png"), width=2000, height=2000, res=300)
          g<-plotPCA(dds, c("batch", "group"), ntop=ntop) + theme_bw3()
          print(g)
          dev.off()
        }
      }

      saveRDS(dds, paste0(outputFilePrefix,curSuffix,".vsd.rds"))

      countNumVsd<-assay(dds)
      colnames(countNumVsd)<-colnames(validCountNum)
      write.table(countNumVsd, paste0(outputFilePrefix,curSuffix,".vsd.txt"),col.names=NA, quote=F, sep="\t")
      
      ylab<-"VSD"
    }
    
    if(!is.na(curgenes)){
      countNumVsd<-countNumVsd[rownames(countNumVsd) %in% curgenes,]
      if("Feature_gene_name" %in% colnames(count)){
        geneNames<-count[rownames(countNumVsd), "Feature_gene_name"]
        if(length(unique(geneNames)) == nrow(countNumVsd)){
          rownames(countNumVsd)<-geneNames
        }
      }
      cat("There are", nrow(countNumVsd), "genes will be used for visualization.\n")
      write.csv(countNumVsd, paste0(outputFilePrefix,curSuffix,".genes.csv"), quote=F)
    }

    hasMultipleGroup<-length(unique(validSampleToGroup$V2)) > 1
    if (hasMultipleGroup) {
      groups<-validSampleToGroup$V2
      if(!draw_all_groups_in_HCA){
        if (is.na(groupColors)){
          unique_groups = unique(groups)
          colors<-makeColors(length(unique_groups))
          names(colors)=unique_groups
        }else{
          colors<-groupColors
        }
        conditionColors<-as.matrix(data.frame(Group=colors[groups]))
      }else{
        cc<-conditionColors[validSampleToGroup$V1,]
        gname = ifelse(title == "all", "Group", title)
        colors<-unique(cc[,gname])
        names(colors)<-unique(validSampleToGroup$V2)
      }
    }else{
      groups<-NA
      colors<-NA
      conditionColors<-NA
    }
    
    #visualization    
    countHT<-countNumVsd

    if(draw_umap){
      library(magrittr)
      library(umap)

      set.seed(20211230)

      normalized_counts<-t(countHT)
      umap_results <- umap::umap(normalized_counts)

      umap_plot_df <- data.frame(umap_results$layout) %>% tibble::rownames_to_column("Sample")
      umap_plot_df$Group<-groups

      g<-ggplot(umap_plot_df, aes(x = X1, y = X2, color = Group)) + geom_point() + theme_bw2() + xlab('UMAP_1') + ylab('UMAP_2')
      for(format in outputFormat){
        if("PDF" == format){
          pdf(paste0(outputFilePrefix,curSuffix,".umap.pdf"),width=7,height=6)
        }else{
          png(paste0(outputFilePrefix,curSuffix,".umap.png"),width=2000,height=1600,res=300)
        }
        print(g)
        dev.off()
      }
    }

    #density plot
    dataForPlot<-reshape2::melt(countHT)
    colnames(dataForPlot)[2]<-"Sample"
    p<-ggplot(dataForPlot, aes(value, colour = Sample)) +
      geom_density() + 
      theme_bw2() + 
      xlab(ylab) + 
      theme(legend.position = "none")

    save_ggplot2_plot(file_prefix=paste0(outputFilePrefix,curSuffix,".density"), 
      outputFormat=outputFormat, 
      width_inch=5, 
      height_inch=3, 
      p)

    nsample=length(unique(dataForPlot$Sample))
    ncol=ceiling(sqrt(nsample))
    nrow=ceiling(nsample/ncol)
    width=max(6, ncol * 2)
    height=max(6, nrow * 2)

    g<-ggplot(dataForPlot) + 
      geom_density(aes(x=value, colour=Sample)) + 
      facet_wrap(~Sample, scales = "free", ncol=ncol) + 
      xlab(ylab) + 
      theme_bw3() + 
      guides(color = FALSE)
    
    save_ggplot2_plot(file_prefix=paste0(outputFilePrefix,curSuffix,".density-individual"), 
      outputFormat=outputFormat, 
      width_inch=width, 
      height_inch=height, 
      g)
   
    if(exists("top25cvInHCA") && top25cvInHCA){
      rv<-rowVars(countNumVsd)
      countVar25<-countNumVsd[rv>=quantile(rv)[4],]
      write.csv(countVar25, paste0(outputFilePrefix,curSuffix,".heatmap.top25variance.csv"))
  
      countList = list(countHT, countVar25)
      names(countList)=c("all", "top25vars")
    }else{
      countList = list(countHT)
      names(countList)=c("all")
    }
    
    cur_name=names(countList)[2]
    for(cur_name in names(countList)){
      cur_counts = countList[[cur_name]]
      
      if (ncol(cur_counts)>1 & nrow(cur_counts)>1) {
        #pca plot
        cat("Drawing PCA for", title, "samples using", nrow(cur_counts), cur_name, "genes.\n")
        gene_suffix = ifelse(cur_name == "all", "", ".top25vars")
        
        drawPCA(paste0(outputFilePrefix, curSuffix, gene_suffix, ".PCA"), cur_counts, showLabelInPCA, groups, colors, outputFormat, width_inch=6, height_inch=5)

        mat_scaled = t(scale(t(cur_counts)))

        if(hasMultipleGroup){
          ha=HeatmapAnnotation( Group=groups,
                                col=list(Group=colors),
                                annotation_legend_param = list(Group = list(ncol = 1, title = "Group", title_position = "topleft")))
        }else{
          ha=NULL
        }

        if(usePearsonInHCA){
          clustering_distance_columns="pearson"
        }else{
          clustering_distance_columns="euclidean"
        }
        ht<-Heatmap( mat_scaled,
                    show_row_names=FALSE,
                    cluster_rows=TRUE, 
                    cluster_columns=TRUE, 
                    show_column_dend=TRUE, 
                    show_row_dend=FALSE,
                    clustering_distance_columns=clustering_distance_columns,
                    name="zscore",
                    top_annotation=ha,
                    use_raster=FALSE)
        file_prefix=paste0(outputFilePrefix, curSuffix, gene_suffix, ".heatmap")
        save_complexheatmap_plot(file_prefix, outputFormat, ht)
      } else {
        cat("Not enough samples or genes. Can't Draw heatmap for", title, "samples.\n")
      }
    }
    
    if (ncol(countNumVsd)>2 & nrow(countNumVsd)>2) {
      cat("Doing correlation analysis of samples ...\n")
      #correlation distribution
      countNumCor<-corTableWithoutZero(countNumVsd,method="spearman")
      write.csv(countNumCor, file=paste0(outputFilePrefix,curSuffix,".Correlation.csv"), row.names=T)

      countNumCorTest<-corTestTableWithoutZero(countNumVsd,method="spearman")
      write.csv(countNumCorTest, file=paste0(outputFilePrefix,curSuffix,".Correlation.Test.csv"), row.names=T)

      if(hasMultipleGroup){
        ha=HeatmapAnnotation( Group=groups,
                              col=list(Group=colors),
                              annotation_legend_param = list(Group = list(ncol = 1, title = "Group", title_position = "topleft")))
      }else{
        ha=NULL
      }

      ht<-Heatmap(countNumCor,
                  show_row_names=TRUE,
                  cluster_rows=TRUE, 
                  cluster_columns=TRUE, 
                  show_column_dend=TRUE, 
                  show_row_dend=FALSE,
                  name="Spearman",
                  top_annotation=ha,
                  use_raster=FALSE)
      file_prefix=paste0(outputFilePrefix, curSuffix, ".Correlation")
      save_complexheatmap_plot(file_prefix, outputFormat, ht)
      
      if (hasMultipleGroup) {
        if(length(unique(validSampleToGroup$V2)) < 3){
          saveInError(paste0("Less than 3 groups. Can't do correlation analysis for group table for ",countTableFile),fileSuffix = paste0(suffix,Sys.Date(),".warning"))
          next
        }
        
        countNumVsdGroup<-mergeTableBySampleGroup(countNumVsd, validSampleToGroup)
        mat_scaled = t(scale(t(countNumVsdGroup)))

        ht<-Heatmap(mat_scaled,
                    show_row_names=FALSE,
                    cluster_rows=TRUE, 
                    cluster_columns=TRUE, 
                    show_column_dend=TRUE, 
                    show_row_dend=FALSE,
                    name="zscore",
                    use_raster=FALSE)
        file_prefix=paste0(outputFilePrefix,curSuffix,".Group.heatmap")
        save_complexheatmap_plot(file_prefix, outputFormat, ht)

        cat("Doing correlation analysis of groups ...\n")
        
        #correlation distribution
        countNumCor<-corTableWithoutZero(countNumVsdGroup,method="spearman")
        write.csv(countNumCor, file=paste0(outputFilePrefix,curSuffix,".Group.Correlation.csv"), row.names=T)
        if (any(is.na(countNumCor)) | length(unique(as.vector(countNumCor)))<=1) {
          cat("NA in group correlation matrix or not enought unique values. Can't draw group correlation figure\n")
          next
        }

        countNumCorTest<-corTestTableWithoutZero(countNumVsdGroup,method="spearman")
        write.csv(countNumCorTest, file=paste0(outputFilePrefix,curSuffix,".Group.Correlation.Test.csv"), row.names=T)
        
        ## Complete cases only, if NA present hclust fails, remove groups with NA
        countNumCor.lower <- countNumCor
        countNumCor.na <- ifelse(lower.tri(countNumCor.lower) ==T, countNumCor, "upper")
        completecases <- rownames(countNumCor)[which(complete.cases(countNumCor.na) ==T)]
        countNumCor <- countNumCor[completecases, completecases]

        ht<-Heatmap(countNumCor,
                    show_row_names=TRUE,
                    cluster_rows=TRUE, 
                    cluster_columns=TRUE, 
                    show_column_dend=TRUE, 
                    show_row_dend=FALSE,
                    name="Spearman",
                    use_raster=FALSE)
        file_prefix=paste0(outputFilePrefix,curSuffix,".Group.Correlation")
        save_complexheatmap_plot(file_prefix, outputFormat, ht)        
      }
    } else {
      cat("Not enough samples or genes. Can't do correlation analysis.\n")
    }
  }
}

if(length(missed_count_tables) == 0){
  writeLines("no count table missing", succeed_file)
  if(file.exists(missed_count_tables_file)){
    unlink(missed_count_tables_file)
  }
}else{
  writeLines(missed_count_tables, "count_table_missing.txt")
}

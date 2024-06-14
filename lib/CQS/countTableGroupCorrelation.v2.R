rm(list=ls()) 
outFile=''
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1=''
parFile2='/nobackup/vickers_lab/projects/20221221_8612_CM_smRNA_human_Glioblastoma/host_genome/bowtie1_genome_1mm_NTA_smallRNA_category/result/CM_8612_Glioblastoma.Category.Table.csv'
parFile3='/nobackup/vickers_lab/projects/20221221_8612_CM_smRNA_human_Glioblastoma/preprocessing/fastqc_post_trim_summary/result/CM_8612_Glioblastoma.countInFastQcVis.Result.Reads.csv'


setwd('/nobackup/vickers_lab/projects/20221221_8612_CM_smRNA_human_Glioblastoma/data_visualization/count_table_correlation_TotalReads/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
# For v2, we use ComplexHeatmap instead of heatmap3, also adjusted the figure size.

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
library(EnhancedVolcano)

suppressPackageStartupMessages(library("ComplexHeatmap"))

myoptions=read_file_map(parSampleFile4, do_unlist=FALSE)

draw_all_groups_in_HCA = is_one(myoptions$draw_all_groups_in_HCA)
draw_umap = is_one(myoptions$draw_umap)
n_first = to_numeric(myoptions$n_first, -1)
usePearsonInHCA = is_one(myoptions$usePearsonInHCA, TRUE)
showLabelInPCA = is_one(myoptions$showLabelInPCA, TRUE)
top25cvInHCA = is_one(myoptions$top25cv_in_hca)
use_green_red_color_in_hca = is_one(myoptions$use_green_red_color_in_hca)
onlySamplesInGroup = is_one(myoptions$onlySamplesInGroup)
hasRowNames = is_one(myoptions$hasRowNames)
useGroupAsBatch = is_one(myoptions$useGroupAsBatch)
transformTable = is_one(myoptions$transformTable)
suffix = to_character(myoptions$suffix, "")
idIndex = to_numeric(myoptions$idIndex, 1)
minMedian = to_numeric(myoptions$minMedian, 1)
minMedianInGroup = to_numeric(myoptions$minMedianInGroup, 1)
fixColorRange = is_one(myoptions$fixColorRange, 1)
useLeastGroups = is_one(myoptions$useLeastGroups, 1)
totalCountKey = to_character(myoptions$totalCountKey, "None")

pca_width_inch=to_numeric(myoptions$pca_width_inch, 4)
pca_height_inch=to_numeric(myoptions$pca_height_inch, 2.5)
pca_point_size=to_numeric(myoptions$pca_point_size, 3)

heatmap_add_width_inch=to_numeric(myoptions$heatmap_add_width_inch, 2)
heatmap_add_height_inch=to_numeric(myoptions$heatmap_add_height_inch, 0)

heatmap_legend_label_fontsize=to_numeric(myoptions$heatmap_legend_label_fontsize, 18)
heatmap_column_name_fontsize=to_numeric(myoptions$heatmap_column_name_fontsize, 18)

legend_label_gp = gpar(fontsize = heatmap_legend_label_fontsize, fontface = "bold")
column_names_gp = gpar(fontsize = heatmap_column_name_fontsize, fontface = "bold")

outputDirectory = to_character(myoptions$outputDirectory, "")
output_include_folder_name = is_one(myoptions$output_include_folder_name, 1)

outputPdf = is_one(myoptions$outputPdf)
outputPng = is_one(myoptions$outputPng, TRUE)
outputTIFF = is_one(myoptions$outputTIFF)

outputFormat<-c()
if(outputPdf){
  outputFormat<-c("PDF")
}
if(outputPng){
  outputFormat<-c(outputFormat, "PNG")
}
if(outputTIFF){
  outputFormat<-c(outputFormat, "TIFF")
}
if(length(outputFormat) == 0){
  outputFormat<-c("PDF")
}

num_top_genes_heatmap=to_numeric(myoptions$num_top_genes_heatmap, 0)

task_suffix<-suffix

countTableFileList<-parSampleFile1
groupFileList<-parSampleFile2
colorFileList<-parSampleFile3

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

error_prefix=paste0(task_suffix, ".", Sys.Date(), ".warning")

#source("/home/zhaos/source/ngsperl/lib/CQS/countTableVisFunctions.R")

##Solving node stack overflow problem start###
#when there are too many genes, drawing dendrogram may failed due to node stack overflow,
#It could be solved by forcing stats:::plotNode to be run as interpreted code rather then byte-compiled code via a nasty hack.
#http://stackoverflow.com/questions/16559250/error-in-heatmap-2-gplots/25877485#25877485

# Convert a byte-compiled function to an interpreted-code function 
unByteCode <- function(fun) {
  FUN <- eval(parse(text=deparse(fun)))
  environment(FUN) <- environment(fun)
  FUN
}

# Replace function definition inside of a locked environment **HACK** 
assignEdgewise <- function(name, env, value) {
  unlockBinding(name, env=env)
  assign( name, envir=env, value=value)
  lockBinding(name, env=env)
  invisible(value)
}

# Replace byte-compiled function in a locked environment with an interpreted-code
# function
unByteCodeAssign <- function(fun) {
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

##Solving node stack overflow problem end###

if (geneFile!="") { #visualization based on genes in geneFile only
  genes<-read.table(geneFile, sep="\t", header=F, stringsAsFactors=F)$V1
  cat("There are", length(genes), "genes in gene file.\n")
}else{
  genes<-NA
}

if(colorFileList != ""){
  groupColors=read_file_map(colorFileList, do_unlist=TRUE)
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

prefix_list=c()

i<-1
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
    next
  }
  
  if (ncol(count)<2) {
    next
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
      if (is.na(groupColors[1])){
        unique_groups = unique(groups)
        colors<-makeColors(length(unique_groups))
        names(colors)=unique_groups
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
    
    cur_file_prefix=paste0(outputFilePrefix, curSuffix)
    prefix_list=c(prefix_list, cur_file_prefix)

    write.table(validSampleToGroup, paste0(cur_file_prefix, ".correlation.groups"), col.names=F, row.names=F, quote=F, sep="\t")
    
    #filter reads/genes by parameter
    validCountNum<-filterCountTable(countNum,validSampleToGroup,minMedian=minMedian,minMedianInGroup=minMedianInGroup)

    if("Feature_gene_biotype" %in% colnames(count)){
      if(length(unique(count$Feature_gene_biotype)) > 1){
        bcounts = count %>% 
          dplyr::select(c("Feature_gene_biotype", colnames(validCountNum))) %>% 
          aggregate(. ~ Feature_gene_biotype, data=., FUN=sum) %>%
          tibble::column_to_rownames("Feature_gene_biotype")
        write.csv(bcounts, paste0(cur_file_prefix, ".biotype_counts.csv"))

        #convert bcounts to percentage table by column
        bperc = t(t(bcounts) / colSums(bcounts) * 100)
        bperc_max = apply(bperc, 1, max)
        cats = bperc_max[bperc_max > 1]
        cats = cats[order(cats, decreasing = T)]
        cat_perc = bperc[names(cats),]
        other_perc = colSums(bperc[!rownames(bperc) %in% names(cats),])
        tperc = rbind(cat_perc, "other"=other_perc)
        write.csv(tperc, paste0(cur_file_prefix, ".biotype_perc.csv"))

        mperc = reshape2::melt(tperc) %>%
          dplyr::rename("Biotype" = "Var1", "Sample"="Var2", "Percentage" = "value")
        mperc$Biotype = factor(mperc$Biotype, levels=rownames(tperc))
        mperc$Sample = factor(mperc$Sample, levels=colnames(validCountNum))

        g=ggplot(mperc, aes(x=Sample, y=Percentage, fill=Biotype)) + 
          geom_bar(stat="identity") + 
          theme_bw3() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                axis.title.x=element_blank())

        width=min(max(7, 0.1 * ncol(validCountNum)) + 3, 50)
        ggsave(paste0(cur_file_prefix, ".biotype_perc.png"), g, width=width, height=6, dpi=300)
      }
    }
    
    #normlize by total reads or VSD
    bNormalizeByCount=FALSE
    if(totalCountFile != ""){
      bNormalizeByCount=totalCountKey != "None"
    }

    if (bNormalizeByCount) { #normlize with total count *10^6
      cat("Normalizing by", totalCountKey, "\n")
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
      write.table(countNumVsd, paste0(cur_file_prefix,".RPM.txt"),col.names=NA, quote=F, sep="\t")
      countNumVsd<-log2(countNumVsd+1)
      ylab<-"log2(Mapped Reads per Million)"
    } else {
      dds=DESeqDataSetFromMatrix(countData = validCountNum, colData = as.data.frame(rep(1,ncol(validCountNum))),design = ~1)
      dds<-try(myEstimateSizeFactors(dds))
      vsdres<-try(dds<-DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
      if (class(vsdres) == "try-error") {
        message=paste0("Warning: varianceStabilizingTransformation function failed.\n",as.character(vsdres))
        warning(message)
        writeLines(message,paste0(cur_file_prefix,".vsd.warning"))
        next;
      }
      
      if(has_batch){
        saveRDS(dds, paste0(cur_file_prefix,".before_removeBatchEffect.dds.rds"))

        countNumVsd<-assay(dds)
        colnames(countNumVsd)<-colnames(validCountNum)
        write.table(countNumVsd, paste0(cur_file_prefix,".before_removeBatchEffect.vsd.txt"),col.names=NA, quote=F, sep="\t")

        dds$batch<-batch_map[colnames(dds)]
        ntop = 3000

        png(paste0(cur_file_prefix,".before_removeBatchEffect.batch.png"), width=2000, height=2000, res=300)
        g<-plotPCA(dds, "batch", ntop=ntop) + theme_bw3()
        print(g)
        dev.off()

        hasMultipleGroup<-length(unique(validSampleToGroup$V2)) > 1
        if(hasMultipleGroup){
          stopifnot(all(rownames(colData(dds)) == validSampleToGroup$V1))
          dds$group = validSampleToGroup$V2

          png(paste0(cur_file_prefix,".before_removeBatchEffect.group.png"), width=2000, height=2000, res=300)
          g<-plotPCA(dds, "group", ntop=ntop) + theme_bw3()
          print(g)
          dev.off()

          png(paste0(cur_file_prefix,".before_removeBatchEffect.batch_group.png"), width=2000, height=2000, res=300)
          g<-plotPCA(dds, c("batch", "group"), ntop=ntop) + theme_bw3()
          print(g)
          dev.off()
        }

        assay(dds) <- limma::removeBatchEffect(assay(dds), dds$batch)
        png(paste0(cur_file_prefix,".after_removeBatchEffect.batch.png"), width=2000, height=2000, res=300)
        g<-plotPCA(dds, "batch", ntop=ntop) + theme_bw3()
        print(g)
        dev.off()

        if(hasMultipleGroup){
          png(paste0(cur_file_prefix,".after_removeBatchEffect.group.png"), width=2000, height=2000, res=300)
          g<-plotPCA(dds, "group", ntop=ntop) + theme_bw3()
          print(g)
          dev.off()

          png(paste0(cur_file_prefix,".after_removeBatchEffect.batch_group.png"), width=2000, height=2000, res=300)
          g<-plotPCA(dds, c("batch", "group"), ntop=ntop) + theme_bw3()
          print(g)
          dev.off()
        }
      }

      saveRDS(dds, paste0(cur_file_prefix,".dds.rds"))

      countNumVsd<-assay(dds)
      colnames(countNumVsd)<-colnames(validCountNum)
      write.table(countNumVsd, paste0(cur_file_prefix,".vsd.txt"),col.names=NA, quote=F, sep="\t")
      
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
      write.csv(countNumVsd, paste0(cur_file_prefix,".genes.csv"), quote=F)
    }

    hasMultipleGroup<-length(unique(validSampleToGroup$V2)) > 1
    if (hasMultipleGroup) {
      groups<-validSampleToGroup$V2
      if(!draw_all_groups_in_HCA){
        if (is.na(groupColors[1])){
          unique_groups = unique(groups)
          colors<-makeColors(length(unique_groups))
          names(colors)=unique_groups
        }else{
          colors<-groupColors
        }
      }else{
        cc<-conditionColors[validSampleToGroup$V1,]
        gname = ifelse(title == "all", "Group", title)
        colors<-unique(cc[,gname])
        names(colors)<-unique(validSampleToGroup$V2)
      }
      conditionColors<-as.matrix(data.frame(Group=colors[groups]))
    }else{
      groups<-NA
      colors<-NA
      conditionColors<-NA
    }

    #density plot
    log2counts<-log2(validCountNum[rownames(countNumVsd), colnames(countNumVsd) ]+1)
    draw_density_plot(log2counts=log2counts, 
                      prefix=paste0(cur_file_prefix,".density"), 
                      outputFormat=outputFormat)

    #visualization    
    if(draw_umap){
      library(magrittr)
      library(umap)

      set.seed(20211230)

      normalized_counts<-t(countNumVsd)
      umap_results <- umap::umap(normalized_counts)

      umap_plot_df <- data.frame(umap_results$layout) %>% tibble::rownames_to_column("Sample")
      umap_plot_df$Group<-groups

      g<-ggplot(umap_plot_df, aes(x = X1, y = X2, color = Group)) + geom_point() + theme_bw2() + xlab('UMAP_1') + ylab('UMAP_2')
      for(format in outputFormat){
        if("PDF" == format){
          pdf(paste0(cur_file_prefix,".umap.pdf"),width=7,height=6)
        }else{
          png(paste0(cur_file_prefix,".umap.png"),width=2000,height=1600,res=300)
        }
        print(g)
        dev.off()
      }
    }

    if(hasMultipleGroup){
      ha=HeatmapAnnotation( Group=groups,
                            col=list(Group=colors),
                            annotation_legend_param = list(Group = list(ncol = 1, 
                                                                        title = "Group", 
                                                                        title_position = "topleft",
                                                                        title_gp=legend_label_gp, 
                                                                        labels_gp=legend_label_gp)))
    }else{
      ha=NULL
    }

    if(num_top_genes_heatmap > 0){
      cur_num = min(num_top_genes_heatmap, nrow(countNumVsd))
      #cur_name=ifelse(bNormalizeByCount, "log(RPM)", "VSD")
      cur_name="Z-Score"

      prefix = paste0(cur_file_prefix, ".top", num_top_genes_heatmap, ".heatmap")

      if(hasMultipleGroup){
        group_levels=levels(groups)
        gname = unique(groups)[1]
        all_genes=c()
        for(gname in group_levels){
          cur_group_vsd<-countNum[, groups == gname]
          cur_group_genes<-rownames(cur_group_vsd)[order(rowMedians(cur_group_vsd),decreasing=T)[1:cur_num]]
          all_genes<-c(all_genes, cur_group_genes)
        }
        all_genes = unique(all_genes)
        cur_vsd<-countNumVsd[all_genes,]
        rownames(cur_vsd)<-gsub(";.+","",rownames(cur_vsd))

        mat_scaled = t(scale(t(cur_vsd)))

        ht_size = draw_heatmap_png( filepath=paste0(prefix, ".png"), 
                                    htdata=mat_scaled, 
                                    name=cur_name, 
                                    save_rds=TRUE,
                                    save_pdf=outputPdf,
                                    add_width_inch=heatmap_add_width_inch,
                                    add_height_inch=heatmap_add_height_inch,
                                    show_row_names=TRUE, 
                                    show_column_names=TRUE,
                                    show_row_dend=TRUE,
                                    column_split = groups,
                                    top_annotation = ha,
                                    column_title=NULL,
                                    column_names_gp = column_names_gp,
                                    legend_gp = legend_label_gp
                                    )
      }else{
        cur_vsd<-countNumVsd[order(rowMedians(countNumVsd),decreasing=T)[1:cur_num],]
        rownames(cur_vsd)<-gsub(";.+","",rownames(cur_vsd))

        mat_scaled = t(scale(t(cur_vsd)))

        ht_size = draw_heatmap_png( filepath=paste0(prefix, ".png"), 
                                    htdata=mat_scaled, 
                                    name=cur_name, 
                                    save_rds=TRUE,
                                    save_pdf=outputPdf,
                                    add_width_inch=heatmap_add_width_inch,
                                    add_height_inch=heatmap_add_height_inch,
                                    show_row_names=TRUE, 
                                    show_column_names=TRUE,
                                    show_row_dend=TRUE,
                                    column_names_gp = column_names_gp,
                                    legend_gp = legend_label_gp
                                    )
      }
    }
   
    if(top25cvInHCA){
      rv<-rowVars(countNumVsd)
      countVar25<-countNumVsd[rv>=quantile(rv)[4],]
      write.csv(countVar25, paste0(cur_file_prefix,".heatmap.top25variance.csv"))
  
      countList = list(countVar25, countNumVsd)
      names(countList)=c("top25vars", "all")
    }else{
      countList = list(countNumVsd)
      names(countList)=c("all")
    }
    
    if(usePearsonInHCA){
      clustering_distance_columns="pearson"
    }else{
      clustering_distance_columns="euclidean"
    }

    cur_name=names(countList)[1]
    for(cur_name in names(countList)){
      cur_counts = countList[[cur_name]]

      show_col_names=ncol(cur_counts) <= 200
      
      if (ncol(cur_counts)>1 & nrow(cur_counts)>1) {
        #pca plot
        cat("Drawing PCA for", title, "samples using", nrow(cur_counts), cur_name, "genes.\n")
        gene_suffix = ifelse(cur_name == "all", "", ".top25vars")
        
        drawPCA(paste0(outputFilePrefix, curSuffix, gene_suffix, ".PCA"), cur_counts, showLabelInPCA, groups, colors, outputFormat, width_inch=pca_width_inch, height_inch=pca_height_inch, point_size=pca_point_size)

        mat_scaled = t(scale(t(cur_counts)))

        draw_heatmap_png( filepath=paste0(outputFilePrefix, curSuffix, gene_suffix, ".heatmap.png"), 
                          htdata=mat_scaled, 
                          name="Z-Score", 
                          save_rds=TRUE,
                          save_pdf=outputPdf,
                          add_width_inch=heatmap_add_width_inch,
                          add_height_inch=heatmap_add_height_inch,
                          show_row_names=FALSE, 
                          show_column_names=show_col_names,
                          show_row_dend=FALSE,
                          top_annotation=ha,
                          clustering_distance_columns=clustering_distance_columns,
                          column_names_gp = column_names_gp,
                          legend_gp = legend_label_gp
                          )
      } else {
        cat("Not enough samples or genes. Can't Draw heatmap for", title, "samples.\n")
      }
    }
    
    if (ncol(countNumVsd)>2 & nrow(countNumVsd)>2) {
      cat("Doing correlation analysis of samples ...\n")
      #correlation distribution
      countNumCor<-corTableWithoutZero(countNumVsd,method="spearman")
      write.csv(countNumCor, file=paste0(cur_file_prefix,".Correlation.csv"), row.names=T)

      countNumCorTest<-corTestTableWithoutZero(countNumVsd,method="spearman")
      write.csv(countNumCorTest, file=paste0(cur_file_prefix,".Correlation.Test.csv"), row.names=T)

      if(any(is.na(countNumCor))){
        saveInError(paste0("NA in correlation matrix. Can't draw correlation figure for ",countTableFile), fileSuffix = error_prefix)
        next
      }
      
      draw_heatmap_png( filepath=paste0(cur_file_prefix, ".Correlation.png"), 
                        htdata=countNumCor, 
                        name="Spearman", 
                        save_rds=TRUE,
                        save_pdf=outputPdf,
                        add_width_inch=heatmap_add_width_inch,
                        add_height_inch=heatmap_add_height_inch,
                        show_row_names=show_col_names, 
                        show_column_names=show_col_names,
                        top_annotation=ha,
                        column_names_gp = column_names_gp,
                        legend_gp = legend_label_gp
                        )
     
      if (hasMultipleGroup) {
        if(length(unique(validSampleToGroup$V2)) < 3){
          saveInError(paste0("Less than 3 groups. Can't do correlation analysis for group table for ", countTableFile),fileSuffix = error_prefix)
          next
        }
        
        countNumVsdGroup<-mergeTableBySampleGroup(countNumVsd, validSampleToGroup)
        mat_scaled = t(scale(t(countNumVsdGroup)))

        draw_heatmap_png( filepath=paste0(cur_file_prefix, ".Group.heatmap.png"), 
                          htdata=mat_scaled, 
                          name="Z-Score", 
                          save_rds=TRUE,
                          save_pdf=outputPdf,
                          add_width_inch=heatmap_add_width_inch,
                          add_height_inch=heatmap_add_height_inch,
                          show_row_names=FALSE, 
                          show_column_names=TRUE,
                          show_row_dend=FALSE,
                          column_names_gp = column_names_gp,
                          legend_gp = legend_label_gp
                          )

        cat("Doing correlation analysis of groups ...\n")
        
        #correlation distribution
        countNumCor<-corTableWithoutZero(countNumVsdGroup,method="spearman")
        write.csv(countNumCor, file=paste0(cur_file_prefix,".Group.Correlation.csv"), row.names=T)
        if (any(is.na(countNumCor)) | length(unique(as.vector(countNumCor)))<=1) {
          cat("NA in group correlation matrix or not enought unique values. Can't draw group correlation figure\n")
          next
        }

        countNumCorTest<-corTestTableWithoutZero(countNumVsdGroup,method="spearman")
        write.csv(countNumCorTest, file=paste0(cur_file_prefix,".Group.Correlation.Test.csv"), row.names=T)
        
        ## Complete cases only, if NA present hclust fails, remove groups with NA
        countNumCor.lower <- countNumCor
        countNumCor.na <- ifelse(lower.tri(countNumCor.lower) ==T, countNumCor, "upper")
        completecases <- rownames(countNumCor)[which(complete.cases(countNumCor.na) ==T)]
        countNumCor <- countNumCor[completecases, completecases]

        draw_heatmap_png( filepath=paste0(cur_file_prefix, ".Group.Correlation.png"), 
                          htdata=countNumCor, 
                          name="Spearman", 
                          save_rds=TRUE,
                          save_pdf=outputPdf,
                          add_width_inch=heatmap_add_width_inch,
                          add_height_inch=heatmap_add_height_inch,
                          show_row_names=TRUE, 
                          show_column_names=TRUE,
                          column_names_gp = column_names_gp,
                          legend_gp = legend_label_gp
                          )
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

writeLines(prefix_list, "prefix_list.txt")

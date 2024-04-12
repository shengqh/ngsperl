rm(list=ls()) 
outFile='9074_ES.host_genome.DESeq2.Matrix'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''
useRawPvalue=1;

setwd('/nobackup/vickers_lab/projects/20230131_9074_ES_ARMseq_human_byMars_hg38/data_visualization/deseq2_host_genome_TotalReads_vis/result')

### Parameter setting end ###

options(bitmapType='cairo')

resultFile<-outFile
deseq2ResultFileList<-parSampleFile1
visLayoutFileList<-parSampleFile2

selectedVars<-c("baseMean","log2FoldChange","pvalue","padj","FoldChange")
#20161207: Don't need this now as significant result file was used to filter significant
#pvalue<-0.05
#foldChange<-1.5

if(!exists("fixColumn")){
  fixColumn<-0
}

#for volcano plot
library(scales)
library(ggplot2)

addVisLayout<-function(datForFigure, visLayoutFileList,LayoutKey="LayoutKey") {
  if (visLayoutFileList!="") {
    visLayout<-read.delim(visLayoutFileList,as.is=T,header=F)
    visLayout<-sapply(split(visLayout[,1],visLayout[,2]),function(x) x)
    if(is.null(ncol(visLayout))){
      visLayout<-data.frame(as.list(visLayout))
    }
    row.names(visLayout)<-visLayout[,"Groups"]
    visLayout<-data.frame(visLayout[,-which(colnames(visLayout)=="Groups")])
    visLayout$Col_Group<-factor(visLayout$Col_Group,levels=unique(visLayout$Col_Group))
    visLayout$Row_Group<-factor(visLayout$Row_Group,levels=unique(visLayout$Row_Group))
    data2Layout<-unique(datForFigure[,LayoutKey])
    for (x in 1:nrow(visLayout)) {
      groupKeys<-strsplit(row.names(visLayout)[x],";")[[1]]
      matchedInd<-1:length(data2Layout)
      for (y in groupKeys) {
        matchedInd<-intersect(matchedInd,grep(y,data2Layout))
      }
      if (length(matchedInd)==1) {
        row.names(visLayout)[x]<-data2Layout[matchedInd]
      } else if (length(matchedInd)>1) {
        message=paste0("Warning: Layout Group: ",row.names(visLayout)[x]," matched to more than one data. The shortest one was used\n Data: \n",paste(data2Layout[matchedInd],collapse="\n"))
        warning(message)
        writeLines(message,paste0(visLayoutFileList,".warning"))
        matchedInd<-matchedInd[which.min(nchar(data2Layout[matchedInd]))]
        row.names(visLayout)[x]<-data2Layout[matchedInd]
      } else {
        message=paste0("Warning: Layout Group: ",row.names(visLayout)[x]," can't match data.\n Data: \n",paste(data2Layout,collapse="\n"))
        writeLines(message,paste0(visLayoutFileList,".warning"))
    warning(message)
      }
    }
    datForFigure<-data.frame(datForFigure,visLayout[datForFigure[,LayoutKey],])
  }
  return(datForFigure)
}

deseq2ResultFile<-read.delim(deseq2ResultFileList,header=F,as.is=T)

deseq2ResultAll<-NULL
for (i in 1:nrow(deseq2ResultFile)) {
  filePath<-deseq2ResultFile[i,1]
  folders<-strsplit(filePath,"\\/")[[1]]
  moduleFolder<-folders[which(folders=="result")-1]
  if(length(moduleFolder) > 1){#multiple result folders in path
    moduleFolder=tail(moduleFolder, 1)
  }
  moduleName<-gsub("_deseq2$","",moduleFolder)
  moduleName<-gsub("^deseq2_","",moduleName)
  moduleName<-gsub("_minicontigs","",moduleName)
  moduleName<-gsub("_contigs","",moduleName)
  if (file.exists(filePath)) {
    deseq2ResultRaw<-read.csv(filePath,header=T,as.is=T)
  } else {
    next;
  }

  deseq2Result<-deseq2ResultRaw[,selectedVars]
  deseq2Result$Module<-moduleName
  deseq2Result$Pairs<-deseq2ResultFile[i,2]
  
  #20161207: Read significant result file so that don't need to filter significant in this module
  fileSigPath<-paste0(tools::file_path_sans_ext(filePath),"_sig.csv")
  deseq2SigResult<-read.csv(fileSigPath,header=T,as.is=T)
  deseq2Result$Significant<-0
  deseq2Result$Significant[which(deseq2ResultRaw[,1] %in% deseq2SigResult[,1])]<-1
  
  deseq2ResultAll<-rbind(deseq2ResultAll,deseq2Result)
}
if(all(grepl("_reads", deseq2ResultAll$Module))){ 
  cat("remove _reads")
  deseq2ResultAll$Module<-gsub("_reads","",deseq2ResultAll$Module)
}
deseq2ResultAll$LayoutKey<-paste0(deseq2ResultAll$Module,"_",deseq2ResultAll$Pairs)

#volcano plot
changeColours<-c(grey="grey",blue="blue",red="red")
diffResult<-as.data.frame(deseq2ResultAll)
diffResult$log10BaseMean<-log10(diffResult$baseMean)
diffResult$colour<-"grey"
#20161207: Use significant result file so that don't need to filter significant in this module
diffResult$colour[which(diffResult$Significant==1 & diffResult$log2FoldChange>0)]<-"red"
diffResult$colour[which(diffResult$Significant==1 & diffResult$log2FoldChange<0)]<-"blue"
#diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange>=log2(foldChange))]<-"red"
#diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange<=-log2(foldChange))]<-"blue"

diffResult$colour=factor(diffResult$colour,levels=c("grey","blue","red"))

diffResult<-addVisLayout(diffResult,visLayoutFileList)

#use log10pvalue instead of pvalue with reverse log scale, otherwise the facet_grid scales might not work.
diffResult$log10pvalue=-log10(diffResult$pvalue)

saveRDS(diffResult, paste0(resultFile, ".rds"))

pair_font_size=16
module_font_size=24

strip_font_family="Times"

get_text_width <- function(txt, font_family, font_size = 10, units = "inches", res=300) {
  tmp_file <- tempfile(fileext = ".png")
  png(tmp_file, res=res)
  par(family = font_family, ps = font_size)
  ret = strwidth(txt, units = units)
  dev.off()
  unlink(tmp_file)

  return(ret)
}

max_pair=max(get_text_width(unique(diffResult$Pairs), font = strip_font_family, font_size = pair_font_size, units = 'inches')) + 0.1
max_module=max(get_text_width(unique(diffResult$Module), font = strip_font_family, font_size = module_font_size, units = 'inches')) + 0.1

width<-max(6,length(unique(diffResult$Pairs))*max_pair)
height<-max(6,length(unique(diffResult$Module))*max_module)

formula=NULL
if((width > height) | fixColumn){
  if (visLayoutFileList!="") {
    formula=as.formula("Row_Group~Col_Group")
  } else {
    formula=as.formula("Module~Pairs")
  }
  xfondsize=pair_font_size
  yfondsize=module_font_size
}else{
  if (visLayoutFileList!="") {
    formula=as.formula("Col_Group~Row_Group")
  } else {
    formula=as.formula("Pairs~Module")
  }
  xfondsize=module_font_size
  yfondsize=pair_font_size
  tmp<-width
  width<-height
  height<-tmp
}

width=width+1
height=height+1.5

diffResult<-diffResult[order(diffResult$colour),]
p<-ggplot(diffResult,aes(x=log2FoldChange,y=log10pvalue))+
    geom_point(aes(colour=colour), size=4)+
    scale_color_manual(values=changeColours,guide="none")+
    scale_x_continuous(name=bquote(log[2](fold~change)),breaks=pretty_breaks(n=4))+
    scale_y_continuous(name=bquote(-log[10](p~value)),breaks=pretty_breaks(n=4))+
    geom_hline(yintercept = 1,colour="grey",linetype = "dotted")+
    geom_vline(xintercept = 0,colour="grey",linetype = "dotted")+
    theme_bw()+
    scale_size(range = c(1, 4))+
    theme(axis.text = element_text(colour = "black",size=20),
        axis.title = element_text(size=20),
        legend.text= element_text(size=20),
        legend.title= element_text(size=20),
        strip.text.x = element_text(family = strip_font_family, size=xfondsize),
        strip.text.y = element_text(family = strip_font_family, size=yfondsize),
        strip.background = element_blank(),
        legend.position="top") +
    facet_grid(formula, scales="free")

ggsave(filename=paste0(resultFile, ".png"), p, width=width, height=height, units="in", dpi=300)
ggsave(filename=paste0(resultFile, ".pdf"), p, width=width, height=height, units="in", dpi=300)

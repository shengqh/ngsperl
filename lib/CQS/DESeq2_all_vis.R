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
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
      log_breaks(base = base), 
      domain = c(1e-100, Inf))
}

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

diffResult<-addVisLayout(diffResult,visLayoutFileList)

if (useRawPvalue==1) {
  p<-ggplot(diffResult,aes(x=log2FoldChange,y=pvalue))+
      scale_y_continuous(trans=reverselog_trans(10),name=bquote(p~value))
} else {
  p<-ggplot(diffResult,aes(x=log2FoldChange,y=padj))+
      scale_y_continuous(trans=reverselog_trans(10),name=bquote(Adjusted~p~value))
}
p<-p+
    geom_point(aes(size=log10BaseMean,colour=colour))+
    scale_color_manual(values=changeColours,guide="none")+
    scale_x_continuous(name=bquote(log[2]~Fold~Change),breaks=pretty_breaks(n=4))+
    geom_hline(yintercept = 1,colour="grey",linetype = "dotted")+
    geom_vline(xintercept = 0,colour="grey",linetype = "dotted")+
    guides(size=guide_legend(title=bquote(log[10]~Base~Mean)))+
    theme_bw()+
    scale_size(range = c(1, 4))+
    theme(axis.text = element_text(colour = "black",size=20),
        axis.title = element_text(size=20),
        legend.text= element_text(size=20),
        legend.title= element_text(size=20))+
    theme(strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        strip.background = element_blank(),
        legend.position="top")
    
width<-max(1600,length(unique(diffResult$Pairs))*800)
height<-max(1600,length(unique(diffResult$Module))*800)

if((width > height) | fixColumn){
  if (visLayoutFileList!="") {
    p<-p+facet_grid(Row_Group~Col_Group,scales = "free")
  } else {
    p<-p+facet_grid(Module~Pairs)
  }
}else{
  if (visLayoutFileList!="") {
    p<-p+facet_grid(Col_Group~Row_Group,scales = "free")
  } else {
    p<-p+facet_grid(Pairs~Module)
  }
  tmp<-width
  width<-height
  height<-tmp
}

png(filename=paste0(resultFile, ".png"), width=width, height=height, res=300)
print(p)
dev.off()


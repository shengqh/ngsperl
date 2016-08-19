resultFile<-outFile
deseq2ResultFileList<-parSampleFile1
visLayoutFileList<-parSampleFile2

selectedVars<-c("baseMean","log2FoldChange","pvalue","padj","FoldChange")
pvalue<-0.05
foldChange<-1.5

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
		row.names(visLayout)<-visLayout[,"Groups"]
		visLayout<-data.frame(visLayout[,-which(colnames(visLayout)=="Groups")])
		visLayout$Col_Group<-factor(visLayout$Col_Group,levels=unique(visLayout$Col_Group))
		visLayout$Row_Group<-factor(visLayout$Row_Group,levels=unique(visLayout$Row_Group))	
		datForFigure<-data.frame(datForFigure,visLayout[datForFigure[,LayoutKey],])
	} else {
		return(datForFigure)
	}
}


deseq2ResultFile<-read.delim(deseq2ResultFileList,header=F,as.is=T)

deseq2ResultAll<-NULL
for (i in 1:nrow(deseq2ResultFile)) {
	filePath<-deseq2ResultFile[i,1]
	folders<-strsplit(filePath,"\\/")[[1]]
	moduleFolder<-folders[which(folders=="result")-1]
	moduleName<-gsub("_deseq2$","",moduleFolder)
	deseq2Result<-read.csv(filePath,header=T,as.is=T)
	deseq2Result<-deseq2Result[,selectedVars]
	deseq2Result$Module<-moduleName
	deseq2Result$Pairs<-deseq2ResultFile[i,2]
	deseq2Result$LayoutKey<-paste0(moduleFolder,"_",deseq2ResultFile[i,2])
	deseq2ResultAll<-rbind(deseq2ResultAll,deseq2Result)
}


#volcano plot
changeColours<-c(grey="grey",blue="blue",red="red")
diffResult<-as.data.frame(deseq2ResultAll)
diffResult$log10BaseMean<-log10(diffResult$baseMean)
diffResult$colour<-"grey"
diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange>=log2(foldChange))]<-"red"
diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange<=-log2(foldChange))]<-"blue"

diffResult<-addVisLayout(diffResult,visLayoutFileList)

p<-ggplot(diffResult,aes(x=log2FoldChange,y=padj))+
		geom_point(aes(size=log10BaseMean,colour=colour))+
		scale_color_manual(values=changeColours,guide = FALSE)+
		scale_y_continuous(trans=reverselog_trans(10),name=bquote(Adjusted~p~value))+
		scale_x_continuous(name=bquote(log[2]~Fold~Change))+
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
				legend.position="top")
		
width<-length(unique(diffResult$Pairs))*800
height<-max(1600,length(unique(diffResult$Module))*800)
png(filename=paste0(resultFile, ".DESeq2.Matrix.png"), width=width, height=height, res=300)
if (visLayoutFileList!="") {
	print(p+facet_grid(Row_Group~Col_Group,scales = "free"))
} else {
	print(p+facet_grid(Module~Pairs))
}
dev.off()






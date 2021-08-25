
###############################################################################
# TODO: Add comment
# 
# Author: zhaos
###############################################################################

deseq2ResultFileTable=parSampleFile1
deseq2ResultFileTable=read.delim(deseq2ResultFileTable,header=F,as.is=T,check.names=F)

###############################
#parameters
###############################
#deseq2ResultTable="/scratch/cqs/zhaos/projects/20171117_BrownRnaSeq/pipeline/deseq2_genetable/result/GroupA_VS_GroupB_DESeq2.csv"
#species='hsa'
#useRawPValue=1
#pCut=0.1

#pCutPathway=pCut
pCutPathway=0.01
pCutPathwayDiff0=pCutPathway/10
useRawPValuePathway=useRawPValue
log2fcColAll=c("log2FoldChange","log2FC","logFC")
if (useRawPValue) {
	pValueColAll=c("pvalue","P.Value")
} else {
	pValueColAll=c("padj","adj.P.Val")
}
if (useRawPValuePathway) {
	pValuePathwayCol="pvalue"
} else {
	pValuePathwayCol="pvalueAdj"
}
if (!exists("geneIdType")) {
  geneIdType="ensembl_gene_id"
}
if (!exists("pathwaysToPlot")) {
  pathwaysToPlot=NULL
}
###############################
#end parameters
###############################



library(KEGGprofile)
library(biomaRt)
library(RCurl)
library(XML)
library(png)
library(ggplot2)
library(KEGG.db)
library(reshape2)

#extract part of color from a color range
col_part<-function(data_all,data_part,col) {
	min_all<-min(data_all,na.rm=T)
	max_all<-max(data_all,na.rm=T)
	min_part<-min(data_part,na.rm=T)
	max_part<-max(data_part,na.rm=T)
	if (min_part<min_all) {
		min_part<-min_all
		warning(paste0("Minimum value in data_part is less than data_all"))
	}
	if (max_part>max_all) {
		max_part<-max_all
		warning(paste0("Maximum value in data_part is larger than data_all"))
	}
	cut_off_low<-round(quantile(1:length(col),(min_part-min_all)/(max_all-min_all)))
	cut_off_high<-round(quantile(1:length(col),(max_part-min_all)/(max_all-min_all)))
	col=col[cut_off_low:cut_off_high]
	return(col)
}

col_by_value<-
  function (x, col, range = NA, breaks = NA, showColorBar = T) 
  {
    if (is.na(range[1])) {
    }
    else {
      x[x < range[1]] <- range[1]
      x[x > range[2]] <- range[2]
    }
    if (is.na(breaks[1])) {
      ff <- seq(min(x,na.rm=TRUE), max(x,na.rm=TRUE), length = length(col))
      bg2 <- apply(as.matrix(as.numeric(unlist(x))), 1, function(x) rank(c(ff, 
                                                                           x), ties.method = "min")[length(col) + 1])
      dens <- matrix(bg2, nrow(x), ncol(x))
      result <- matrix(col[dens], nrow = nrow(x), ncol = ncol(x))
      row.names(result) <- row.names(x)
      if (showColorBar) {
        image(x = 1, y = as.matrix(ff), z = t(ff), col = col, 
              xaxt = "n", ylab = "")
        box()
      }
      return(result)
    }
    else {
      temp <- cut(as.numeric(unlist(x)), breaks = breaks, include.lowest = T)
      if (length(col) != length(levels(temp))) {
        stop("length:col != length: cut result")
      }
      result <- matrix(col[as.numeric(temp)], nrow = nrow(x), 
                       ncol = ncol(x))
      row.names(result) <- row.names(x)
      if (showColorBar) {
        par(mar = c(5, 9, 2, 2))
        image(x = 1, y = as.matrix(1:(length(breaks) - 1)), 
              z = t(1:(length(breaks) - 1)), col = col, xaxt = "n", 
              yaxt = "n", ylab = "", xlab = "")
        axis(2, at = 1:(length(breaks) - 1), labels = levels(temp), 
             las = 1)
        box()
      }
      return(result)
    }
  }

plot_pathway_overall<-
function (gene_values, species = "hsa", pathwayNumInFigure = 5, 
    rankByVar = colnames(gene_values)[1]) 
{
    kegg_enriched_pathway = find_enriched_pathway(row.names(gene_values), 
        species = species, returned_pvalue = 1, returned_adjpvalue = 1, 
        returned_genenumber = 5)
    names(kegg_enriched_pathway[[2]]) = make.names(kegg_enriched_pathway[[1]]$Pathway_Name)
    geneValuesInPathway <- lapply(kegg_enriched_pathway[[2]], 
        function(x) gene_values[intersect(x, row.names(gene_values)), 
            , drop = FALSE])
    dataForPlot = NULL
    for (i in 1:length(geneValuesInPathway)) {
        dataForPlot = rbind(dataForPlot, cbind(geneValuesInPathway[[i]], 
            Pathway = names(geneValuesInPathway)[i],Id=row.names(kegg_enriched_pathway[[1]])[i]))
    }
    dataForPlot = reshape2::melt(dataForPlot, id.vars = c("Pathway","Id"), 
        variable.name = "Sample", value.name = "value")
    temp = dataForPlot[which(dataForPlot[, "Sample"] == rankByVar), 
        ]
    pathwayOrder = tapply(temp$value, temp$Pathway, median, na.rm = T)
    pathwayOrder = names(pathwayOrder)[order(pathwayOrder)]
    if (length(pathwayOrder) > pathwayNumInFigure * 2) {
        pathwayOrder = unique(c(head(pathwayOrder, pathwayNumInFigure), 
            tail(pathwayOrder, pathwayNumInFigure)))
    }
    dataForPlot = dataForPlot[which(dataForPlot$Pathway %in% 
        pathwayOrder), ]
    dataForPlot$Pathway = factor(dataForPlot$Pathway, levels = pathwayOrder)
    p = ggplot(dataForPlot, aes_string(x = "Pathway", y = "value")) + 
        geom_boxplot(aes_string(colour = "Sample"))
    return(p + coord_flip() + theme(axis.text.x = element_text(angle = 90, 
        hjust = 1)))
}
###############################
#end functions
###############################


###############################
#figure for all comparisions
###############################
resultTable=NULL
resultTableList=list()
for (j in 1:nrow(deseq2ResultFileTable)) {
  deseq2ResultTable=deseq2ResultFileTable[j,1]
  resultTableOne<-read.csv(deseq2ResultTable,header=T,as.is=T,row.names=1)
  row.names(resultTableOne)=gsub("\\.\\d+$","",row.names(resultTableOne))
  resultTableOne$rownames=row.names(resultTableOne)
  resultTableOne$Comparison=tools::file_path_sans_ext(basename(deseq2ResultTable))
  resultTableList[[j]]=resultTableOne
#  resultTable=rbind(resultTable,resultTableOne)
}
colToMerge=lapply(resultTableList,colnames)
colToMerge=table(unlist(colToMerge))
colToMerge=names(colToMerge[which(colToMerge==max(colToMerge))])
for (j in 1:length(resultTableList)) {
  resultTable=rbind(resultTable,resultTableList[[j]][,colToMerge])
}

log2fcCol=intersect(colnames(resultTable),log2fcColAll)
pValueCol=intersect(colnames(resultTable),pValueColAll)
if (length(log2fcCol)!=1 | length(pValueCol)!=1) {
  warning(paste0("Can't find log2fcCol or pValueCol. Please check colnames in ",deseq2ResultTable));
  next;
}

#change gene expression
#To remove outlier. For not significant genes, change their fold change to less than 2.
temp<-resultTable[,c(log2fcCol,"rownames","Comparison"),drop=FALSE]
for (i in which(resultTable[,pValueCol]>pCut)) {
  if (temp[i,1]>=0) {
    temp[i,1]<-min(log2(2),temp[i,1])
  } else {
    temp[i,1]<-max(-log2(2),temp[i,1])
  }
}
temp=reshape2::dcast(temp,rownames~Comparison,value.var=log2fcCol)
row.names(temp)=temp$rownames
temp=temp[,-1,drop=FALSE]
#	head(temp)
if (species=="hsa") {
  dataset="hsapiens_gene_ensembl"
} else if (species=="mmu") {
  dataset="mmusculus_gene_ensembl"
} else if (species=="rno") {
  dataset="rnorvegicus_gene_ensembl"
} else {
  stop(paste0("species only supports hsa, mmu, or rno at this time."))
}
resultTableFcToGene<-convertId(temp,filters=geneIdType,dataset=dataset)
#	geneExpr<-resultTableFcToGene[,1]
#	names(geneExpr)<-row.names(resultTableFcToGene)
png(paste0(outFile,".KEGG.OverallExpression.png"),width=3000,height=4000,res=300)
p=plot_pathway_overall(resultTableFcToGene,species = species,pathwayNumInFigure = 10)
plot(p + theme(legend.position="top"))
dev.off()
#keggEnrichedPathway<-find_enriched_pathway(row.names(resultTableFcToGene),species=species,returned_genenumber = 5,returned_pvalue=1,returned_adjpvalue = 1)
#dim(KEGGresult1[[1]])
pathwaysToPlot=c(pathwaysToPlot,unique(p$data$Id))

###################################################
#If any significant pathway, plot gene expression changes in it
###################################################	
resultTablePSig=resultTable[which(resultTable[,pValueCol]<=pCut),pValueCol,drop=FALSE]
if (nrow(resultTablePSig)>10) { #at least 10 genes to do enrichment, or meaningless
  temp=1-resultTablePSig #1-p value, so that we can keep the gene with smallest p value when using genesKept="abs"
  resultTablePSigToGene<-convertId(temp,filters=geneIdType,dataset=dataset,genesKept="abs",keepNoId=FALSE)
  genesEnz<-na.omit(row.names(resultTablePSigToGene))
  keggSigGeneEnrichedPathway<-find_enriched_pathway(genesEnz,species=species,returned_genenumber = 1,returned_pvalue=1,returned_adjpvalue = 1)
  
} else {
  warning(paste0("Less than 10 significant genes in ",deseq2ResultTable," Skip.")) 
  next;
}

#differential genes enriched pathways
sigPathwayInd=which(keggSigGeneEnrichedPathway[[1]][,pValuePathwayCol]<=pCutPathway)
if (length(sigPathwayInd)>0) {
  pathwaysToPlot=unique(c(pathwaysToPlot,row.names(keggSigGeneEnrichedPathway[[1]])[sigPathwayInd]))
}

if (length(pathwaysToPlot)>0) {
  #In case of unbalance fold change. balance the color
  colPart<-col_part(c(-2.5,2.5), range(resultTableFcToGene[,1:ncol(resultTableFcToGene)],na.rm=T),col=colorRampPalette(c('green','black','red'))(1024))
  png(paste0(outFile,"_","FoldChange_colorBar",".png"),width=300)
  col<-col_by_value(resultTableFcToGene,col=colPart,range=c(-3,3))
  dev.off()
  
  #make border color to show significant genes
  resultTableFcToGeneBorderCol=matrix("grey",ncol=ncol(resultTableFcToGene),nrow=nrow(resultTableFcToGene))
  row.names(resultTableFcToGeneBorderCol)=row.names(resultTableFcToGene)
  resultTableFcToGeneBorderCol[row.names(resultTablePSigToGene),]="darkorange1"
  
  for (pathway_id in pathwaysToPlot) {
    #			pathway_id=i
    #			download_KEGGfile(pathway_id)
    result_name=paste0(outFile,"_",species,"_",pathway_id,"_fc",".png")
    temp<-try(plot_pathway(resultTableFcToGene,type="bg",bg_col=col,text_col="white",magnify=1.4,species=species,pathway_id=pathway_id,genes_kept="abs",result_name=result_name,border_col=resultTableFcToGeneBorderCol))
  }
} else {
  warning(paste0("No pathway has significant gene enrichment! Can't make pathway gene expression figures."))
}


########################
#generate report
########################
#make files to include all keggprofile result, to be used in the Rmd file which is designed for limma pipeline
keggprofilePathwayFigures=list.files(getwd(),pattern = "_fc.png$",full.names = TRUE)
keggprofileOverallFigures=list.files(getwd(),pattern = ".KEGG.OverallExpression.png$",full.names = TRUE)
files=cbind(c(keggprofilePathwayFigures,keggprofileOverallFigures))
#render report
keggReportFile=paste0(basename(keggprofileOverallFigures[1]),".Report.html")
rmarkdown::render("KEGGprofileReport.Rmd",output_file=keggReportFile,output_dir=getwd())


# 
# 
# ###############################
# #figure for each comparision
# ###############################
# for (j in 1:nrow(deseq2ResultFileTable)) {
# 	deseq2ResultTable=deseq2ResultFileTable[j,1]
# 	keggOutFileName<-paste0(tools::file_path_sans_ext(basename(deseq2ResultTable)),"_KEGG")
# 	dir.create(keggOutFileName,showWarnings = FALSE)
# 	
# 	resultTable<-read.csv(deseq2ResultTable,header=T,as.is=T,row.names=1)
# 	row.names(resultTable)=gsub("\\.\\d+$","",row.names(resultTable))
# 	log2fcCol=intersect(colnames(resultTable),log2fcColAll)
# 	pValueCol=intersect(colnames(resultTable),pValueColAll)
# 	if (length(log2fcCol)!=1 | length(pValueCol)!=1) {
# 	  warning(paste0("Can't find log2fcCol or pValueCol. Please check colnames in ",deseq2ResultTable));
# 	  next;
# 	}
# 	
# 	
# #change gene expression
# #To remove outlier. For not significant genes, change their fold change to less than 2.
# 	temp<-resultTable[,c(log2fcCol),drop=FALSE]
# 	for (i in which(resultTable[,pValueCol]>pCut)) {
# 		if (temp[i,1]>=0) {
# 			temp[i,1]<-min(log2(2),temp[i,1])
# 		} else {
# 			temp[i,1]<-max(-log2(2),temp[i,1])
# 		}
# 	}
# #	head(temp)
# 	if (species=="hsa") {
# 		dataset="hsapiens_gene_ensembl"
# 	} else if (species=="mmu") {
# 		dataset="mmusculus_gene_ensembl"
# 	} else if (species=="rno") {
# 		dataset="rnorvegicus_gene_ensembl"
# 	} else {
# 		stop(paste0("species only supports hsa, mmu, or rno at this time."))
# 	}
# 	resultTableFcToGene<-convertId(temp,filters=geneIdType,dataset=dataset)
# #	geneExpr<-resultTableFcToGene[,1]
# #	names(geneExpr)<-row.names(resultTableFcToGene)
# 	png(paste0(keggOutFileName,".OverallExpression.png"),width=1500,height=3000,res=300)
# 	plot_pathway_overall(resultTableFcToGene,species = species)
# 	dev.off()
# 	keggEnrichedPathway<-find_enriched_pathway(row.names(resultTableFcToGene),species=species,returned_genenumber = 5,returned_pvalue=1,returned_adjpvalue = 1)
# 	#dim(KEGGresult1[[1]])
# 	
# 	
# 	###################################################
# 	#If any significant pathway, plot gene expression changes in it
# 	###################################################	
# 	resultTablePSig=resultTable[which(resultTable[,pValueCol]<=pCut),pValueCol,drop=FALSE]
# 	if (nrow(resultTablePSig)>10) { #at least 10 genes to do enrichment, or meaningless
# 	  temp=1-resultTablePSig #1-p value, so that we can keep the gene with smallest p value when using genesKept="abs"
# 	  resultTablePSigToGene<-convertId(temp,filters=geneIdType,dataset=dataset,genesKept="abs",keepNoId=FALSE)
# 	  genesEnz<-na.omit(row.names(resultTablePSigToGene))
# 	  keggSigGeneEnrichedPathway<-find_enriched_pathway(genesEnz,species=species,returned_genenumber = 1,returned_pvalue=1,returned_adjpvalue = 1)
# 	  
# 	} else {
# 	   warning(paste0("Less than 10 significant genes in ",deseq2ResultTable," Skip.")) 
# 	   next;
# 	}
# 
# 	#differential genes enriched pathways
# 	sigPathwayInd=which(keggSigGeneEnrichedPathway[[1]][,pValuePathwayCol]<=pCutPathway)
# 	if (length(sigPathwayInd)>0) {
# 	  #In case of unbalance fold change. balance the color
# 	  colPart<-col_part(c(-2.5,2.5), range(resultTableFcToGene[,1]),col=colorRampPalette(c('green','black','red'))(1024))
# 	  png(paste0(keggOutFileName,"/","FoldChange_colorBar",".png"),width=300)
# 	  col<-col_by_value(resultTableFcToGene,col=colPart,range=c(-3,3))
# 	  dev.off()
# 	  
# 	  #make border color to show significant genes
# 	  resultTableFcToGeneBorderCol=matrix("grey",ncol=1,nrow=nrow(resultTableFcToGene))
# 	  row.names(resultTableFcToGeneBorderCol)=row.names(resultTableFcToGene)
# 	  resultTableFcToGeneBorderCol[row.names(resultTablePSigToGene),1]="darkorange1"
# 	  
# 	  for (pathway_id in row.names(keggSigGeneEnrichedPathway[[1]])[sigPathwayInd]) {
# 	    #			pathway_id=i
# 	    #			download_KEGGfile(pathway_id)
# 	    result_name=paste0(keggOutFileName,"/",species,"_",pathway_id,"_fc",".png")
# 	    temp<-try(plot_pathway(resultTableFcToGene,type="bg",bg_col=col,text_col="white",magnify=1.4,species=species,pathway_id=pathway_id,genes_kept="abs",result_name=result_name,border_col=resultTableFcToGeneBorderCol))
# 	  }
# 	} else {
# 	  warning(paste0("No pathway has significant gene enrichment! Can't make pathway gene expression figures."))
# 	}
# 	
# 	
# 	##########################################
# #output pathway information table
# 	##########################################
# 
# 	keggResultOut<-data.frame(Pathway_Id=row.names(keggEnrichedPathway[[1]]),stringsAsFactors=FALSE)
# 	keggResultOut$Pathway_Name<-keggEnrichedPathway[[1]][,c("Pathway_Name")]
# 	keggResultOut$Pathway_SignificantDiffGene<-keggSigGeneEnrichedPathway[[1]][keggResultOut$Pathway_Id,c("Gene_Found")]
# 	keggResultOut$Pathway_DataGene<-keggEnrichedPathway[[1]][,c("Gene_Found")]
# 	keggResultOut$Pathway_AllGene<-keggEnrichedPathway[[1]][,c("Gene_Pathway")]
# #	keggResultOut$Pathway_Log2FoldChange<-pathwayMedian[keggResultOut$Pathway_Name]
# #	keggResultOut$Pathway_Log2FoldChangePValue<-pathwayDiff0[keggResultOut$Pathway_Name]
# 	keggResultOut$Pathway_SignificantDiffGenePValue<-keggSigGeneEnrichedPathway[[1]][keggResultOut$Pathway_Id,c("pvalue")]
# 	keggResultOut$Pathway_SignificantDiffGeneAdjPValue<-keggSigGeneEnrichedPathway[[1]][keggResultOut$Pathway_Id,c("pvalueAdj")]
# 	
# 	keggTableFileName<-paste0(keggOutFileName,".csv")
# 	write.csv(keggResultOut,keggTableFileName,row.names=FALSE)
# 	
# }
# 

save.image(paste0(outFile,".KEGGprofile.RData"))




###############################################################################
# TODO: Add comment
# 
# Author: zhaos
###############################################################################

deseq2ResultFileTable=parSampleFile1
deseq2ResultFileTable=read.delim(deseq2ResultFileTable,header=F,as.is=T)

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
if (useRawPValue) {
	pValueCol="pvalue"
} else {
	pValueCol="padj"
}
if (useRawPValuePathway) {
	pValuePathwayCol="pvalue"
} else {
	pValuePathwayCol="pvalueAdj"
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

###############################
#end functions
###############################



for (j in 1:nrow(deseq2ResultFileTable)) {
	deseq2ResultTable=deseq2ResultFileTable[j,1]
	keggOutFileName<-paste0(tools::file_path_sans_ext(basename(deseq2ResultTable)),"_KEGG")
	dir.create(keggOutFileName,showWarnings = FALSE)
	
	resultTable<-read.csv(deseq2ResultTable,header=T,as.is=T,row.names=1)
	row.names(resultTable)=gsub("\\.\\d+$","",row.names(resultTable))

#change gene expression
#To remove outlier. For not significant genes, change their fold change to less than 2.
	temp<-resultTable[,c("log2FoldChange"),drop=FALSE]
	for (i in which(resultTable[,pValueCol]>pCut)) {
		if (temp[i,1]>=0) {
			temp[i,1]<-min(log2(2),temp[i,1])
		} else {
			temp[i,1]<-max(-log2(2),temp[i,1])
		}
	}
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
	resultTableFcToGene<-convertId(temp,filters="ensembl_gene_id",dataset=dataset)
#	geneExpr<-resultTableFcToGene[,1]
#	names(geneExpr)<-row.names(resultTableFcToGene)
	png(paste0(keggOutFileName,".OverallExpression.png"),width=1500,height=3000,res=300)
	plot_pathway_overall(resultTableFcToGene,species = species)
	dev.off()
	keggEnrichedPathway<-find_enriched_pathway(row.names(resultTableFcToGene),species=species,returned_genenumber = 5,returned_pvalue=1,returned_adjpvalue = 1)
	#dim(KEGGresult1[[1]])
	
	
	###################################################
	#If any significant pathway, plot gene expression changes in it
	###################################################	
	resultTablePSig=resultTable[which(resultTable[,pValueCol]<=pCut),pValueCol,drop=FALSE]
	if (nrow(resultTablePSig)>10) { #at least 10 genes to do enrichment, or meaningless
	  temp=1-resultTablePSig #1-p value, so that we can keep the gene with smallest p value when using genesKept="abs"
	  resultTablePSigToGene<-convertId(temp,filters="ensembl_gene_id",dataset=dataset,genesKept="abs",keepNoId=FALSE)
	  genesEnz<-na.omit(row.names(resultTablePSigToGene))
	  keggSigGeneEnrichedPathway<-find_enriched_pathway(genesEnz,species=species,returned_genenumber = 1,returned_pvalue=1,returned_adjpvalue = 1)
	  
	} else {
	   warning(paste0("Less than 10 significant genes in ",deseq2ResultTable," Skip.")) 
	   next;
	}

	#differential genes enriched pathways
	sigPathwayInd=which(keggSigGeneEnrichedPathway[[1]][,pValuePathwayCol]<=pCutPathway)
	if (length(sigPathwayInd)>0) {
	  #In case of unbalance fold change. balance the color
	  colPart<-col_part(c(-2.5,2.5), range(resultTableFcToGene[,1]),col=colorRampPalette(c('green','black','red'))(1024))
	  png(paste0(keggOutFileName,"/","FoldChange_colorBar",".png"),width=300)
	  col<-col_by_value(resultTableFcToGene,col=colPart,range=c(-3,3))
	  dev.off()
	  
	  #make border color to show significant genes
	  resultTableFcToGeneBorderCol=matrix("grey",ncol=1,nrow=nrow(resultTableFcToGene))
	  row.names(resultTableFcToGeneBorderCol)=row.names(resultTableFcToGene)
	  resultTableFcToGeneBorderCol[row.names(resultTablePSigToGene),1]="darkorange1"
	  
	  for (pathway_id in row.names(keggSigGeneEnrichedPathway[[1]])[sigPathwayInd]) {
	    #			pathway_id=i
	    #			download_KEGGfile(pathway_id)
	    result_name=paste0(keggOutFileName,"/",species,"_",pathway_id,"_fc",".png")
	    temp<-plot_pathway(resultTableFcToGene,type="bg",bg_col=col,text_col="white",magnify=1.4,species=species,pathway_id=pathway_id,genes_kept="abs",result_name=result_name,border_col=resultTableFcToGeneBorderCol)
	  }
	} else {
	  warning(paste0("No pathway has significant gene enrichment! Can't make pathway gene expression figures."))
	}
	
	
	##########################################
#output pathway information table
	##########################################

	keggResultOut<-data.frame(Pathway_Id=row.names(keggEnrichedPathway[[1]]),stringsAsFactors=FALSE)
	keggResultOut$Pathway_Name<-keggEnrichedPathway[[1]][,c("Pathway_Name")]
	keggResultOut$Pathway_SignificantDiffGene<-keggSigGeneEnrichedPathway[[1]][keggResultOut$Pathway_Id,c("Gene_Found")]
	keggResultOut$Pathway_DataGene<-keggEnrichedPathway[[1]][,c("Gene_Found")]
	keggResultOut$Pathway_AllGene<-keggEnrichedPathway[[1]][,c("Gene_Pathway")]
#	keggResultOut$Pathway_Log2FoldChange<-pathwayMedian[keggResultOut$Pathway_Name]
#	keggResultOut$Pathway_Log2FoldChangePValue<-pathwayDiff0[keggResultOut$Pathway_Name]
	keggResultOut$Pathway_SignificantDiffGenePValue<-keggSigGeneEnrichedPathway[[1]][keggResultOut$Pathway_Id,c("pvalue")]
	keggResultOut$Pathway_SignificantDiffGeneAdjPValue<-keggSigGeneEnrichedPathway[[1]][keggResultOut$Pathway_Id,c("pvalueAdj")]
	
	keggTableFileName<-paste0(keggOutFileName,".csv")
	write.csv(keggResultOut,keggTableFileName,row.names=FALSE)
	
}


save.image(paste0(outFile,".KEGGprofile.RData"))



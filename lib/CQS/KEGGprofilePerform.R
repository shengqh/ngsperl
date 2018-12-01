# TODO: Add comment
# 
# Author: zhaos
###############################################################################

#setwd("/scratch/cqs/zhaos/temp")

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

###############################
#functions
###############################
#source("D:\\source\\r_cqs\\myPkg\\R\\convertIdOneToOne.R")
#source("D:\\source\\r_cqs\\myPkg\\R\\convertId.R")
#source("D:\\source\\r_cqs\\myPkg\\R\\colorPart.R")
convertId<-function(x,dataset="hsapiens_gene_ensembl",filters="uniprotswissprot",attributes =c(filters,"entrezgene"),genesKept=c('foldchange','first','random','var','abs'),keepNoId=T,keepMultipleId=F,verbose=F) {
#	if (! require("biomaRt")) {
#		cat("biomaRt package is needed but not installed in this computer. Will install it from bioconductor.\n")
#		flush.console()
#		source("http://bioconductor.org/biocLite.R")
#		biocLite("biomaRt")
#		if (!require(biomaRt)) {stop("Package biomaRt can't be installed")}
#	}
	if (missing(genesKept)) {
		genesKept<-"var"
	} else {
		genesKept<-match.arg(genesKept)
	}
	if (verbose) {
		cat("Now conectting with ensembl. Internet acess is needed and it may use 30 seconds.\n")
		flush.console()
	}
	
	#temp
	#ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
	#temp<-listAttributes(ensembl)
	ensembl = useMart("ensembl",dataset=dataset)
	
	oldId<-row.names(x)
	newIdTable<-getBM(attributes =attributes,filters=filters,values=oldId,mart = ensembl)
	newIdTable<-newIdTable[which(newIdTable[,1]!="" & newIdTable[,2]!=""),]
	
	temp1<-which(oldId %in% newIdTable[,1])
	temp2<-nrow(x)-length(temp1)
	xNoId<-NULL
	if (keepNoId) {
		if (verbose) {
			cat(paste("No ID, Keep: ",temp2," genes can't find their ",attributes[2]," ID. They will be attched at the end of data with their original ID.\n",sep=""))
		}
		xNoId<-x[-temp1,,drop=F]
	} else {
		if (verbose) {
			cat(paste("No ID, Discard: ",temp2," genes can't find their ",attributes[2]," ID. They will be discard.\n",sep=""))
		}
	}
	
	temp2<-split(newIdTable,newIdTable[,1])
	if (keepMultipleId) {
		newIdTable<-sapply(temp2,function(x) return(paste(x[,2],collapse=";")))
		if (verbose) {
			cat(paste("Multiple IDs, Keep: ",length(which(sapply(temp2,function(x) nrow(x))>=2))," genes have more than one ",attributes[2]," IDs. All of these ",attributes[2]," IDs will be stored.\n",sep=""))
		}
	} else {
		newIdTable<-sapply(temp2,function(x) return(x[1,2]))
		if (verbose) {
			cat(paste("Multiple IDs, Discard: ",length(which(sapply(temp2,function(x) nrow(x))>=2))," genes have more than one ",attributes[2]," IDs. Only the first ",attributes[2]," ID for each gene will be used.\n",sep=""))
		}
	}
	names(newIdTable)<-names(temp2)
	
	result<-newIdMatrix(x,genesKept=genesKept,convertIdTable=newIdTable)
	
	if (keepMultipleId) {
		temp<-strsplit(row.names(result),";")
		temp1<-sapply(temp,length)
		temp2<-unlist(temp)
		result<-result[rep(1:length(temp1),temp1),]
		row.names(result)<-temp2
	}
	if (keepNoId) {
		result<-rbind(result,xNoId)
	}
	return(result)
}

newIdMatrix<-function(x,convertIdTable,genesKept=c("var","foldchange","abs","first","random")) {
	convertIdTable<-convertIdTable[which(convertIdTable!="" & names(convertIdTable)!="")]
	x<-x[names(convertIdTable),,drop=F]
	if (missing(genesKept)) {
		genesKept<-"var"
	} else {
		genesKept<-match.arg(genesKept)
	}
	if (genesKept=="foldchange") {
		temp<-apply(x,1,range,na.rm=T)
		testStat<-temp[2,]-temp[1,]
	} else if (genesKept=="first") {
		testStat<-rep(1,nrow(x))
		names(testStat)<-row.names(x)
	} else if (genesKept=="random") {
		testStat<-sample(1:nrow(x),nrow(x))
		names(testStat)<-row.names(x)
	} else if (genesKept=="var") {
		testStat<-apply(x,1,var,na.rm=T)
	} else if (genesKept=="abs") {
		testStat<-apply(x,1,function(y) max(abs(y),na.rm=TRUE))
	}
	testStat[is.na(testStat)]<--Inf #Some temp has a NA value if genesKept=="var" and only one sample has value
	temp<-split(testStat,convertIdTable)
	result<-x[unlist(sapply(temp, function(y) names(which.max(y)))),,drop=FALSE]
	row.names(result)<-names(temp)
	return(result)
}

convertIdOneToOne<-function(x,dataset="hsapiens_gene_ensembl",filters="uniprotswissprot",attributes =c(filters,"entrezgene"),verbose=FALSE) {
	if (verbose) {
		cat("Now conectting with ensembl. Internet acess is needed and it may use 30 seconds.\n")
		flush.console()
	}
	
	ensembl = useMart("ensembl",dataset=dataset)
	newIdTable<-getBM(attributes =attributes,filters=filters,values=x,mart = ensembl)
	newIdTable<-newIdTable[which(newIdTable[,1]!="" & newIdTable[,2]!=""),]
	if (verbose & any(table(newIdTable[,1])>1)) {
		cat(paste0(length(which(table(newIdTable[,1])>1))," Ids can be converted to more than one ",attributes[2],", will use the first one.\n"))
	}
	result<-newIdTable[match(x, newIdTable[,1]),2]
	if (verbose & any(is.na(result))) {
		cat(paste0(length(which(is.na(result)))," Ids cant't be converted to ",attributes[2],", will be set as NA.\n"))
	}
	names(result)<-x
	return(result)
}

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
#To remove outlier. For not significant genes, change their fold change to less than 1.5.
	temp<-resultTable[,c("log2FoldChange"),drop=FALSE]
	for (i in which(resultTable[,pValueCol]>pCut)) {
		if (temp[i,1]>=0) {
			temp[i,1]<-min(log2(1.5),temp[i,1])
		} else {
			temp[i,1]<-max(log2(2/3),temp[i,1])
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
	geneExpr<-resultTableFcToGene[,1]
	names(geneExpr)<-row.names(resultTableFcToGene)
	keggEnrichedPathway<-find_enriched_pathway(names(geneExpr),species=species,returned_genenumber = 5,returned_pvalue=1,returned_adjpvalue = 1)
#dim(KEGGresult1[[1]])
	
	genes<-row.names(resultTable)[which(resultTable[,pValueCol]<=pCut)]
	genesEnz<-convertIdOneToOne(genes,filters="ensembl_gene_id",dataset=dataset)
	genesEnz<-na.omit(genesEnz)
	keggSigGeneEnrichedPathway<-find_enriched_pathway(genesEnz,species=species,returned_genenumber = 1,returned_pvalue=1,returned_adjpvalue = 1)
	
	
	
	###################################################
#expression changes for all genes in each pathway
	###################################################
#make gene fold changes for each pathway
	dataForPlot<-data.frame(stringsAsFactors=FALSE)
	for (i in 1:nrow(keggEnrichedPathway[[1]])) {
		tempGeneExp<-geneExpr[keggEnrichedPathway[[2]][[i]]]
		tempPathwayName<-keggEnrichedPathway[[1]][i,1]
		dataForPlot<-rbind(dataForPlot,cbind(tempPathwayName,tempGeneExp))
	}
#head(dataForPlot)
#dim(dataForPlot)
	colnames(dataForPlot)<-c("Pathway","FoldChange")
	dataForPlot[,2]<-as.numeric(as.character(dataForPlot[,2]))
	
#plot boxplot for pathway fold changes. Only show pathways with fold changes different than 0 (t test p<=0.01)
	pathwayMedian<-tapply(dataForPlot[,2],dataForPlot[,1],median)
	pathwayMaxChange<-tapply(dataForPlot[,2],dataForPlot[,1],function(x) max(abs(x)))
	pathwayDiff0<-tapply(dataForPlot[,2],dataForPlot[,1],function(x) t.test(x,mu=0)$p.value)
	selectedPathways<-names(which(pathwayDiff0<=pCutPathwayDiff0))
	if (length(selectedPathways)>0) {
		selectedPathways<-selectedPathways[which(selectedPathways!="Cyanoamino acid metabolism")] #not in human
		selectedPathways<-selectedPathways[which(selectedPathways!="Metabolic pathways")] #not in human
		dataForPlot1<-dataForPlot[which(dataForPlot[,1] %in% selectedPathways),]
		pathwayMedian1<-pathwayMedian[selectedPathways]
		dataForPlot1[,1]<-factor(as.character(dataForPlot1[,1]),levels=names(pathwayMedian1)[order(pathwayMedian1)])
		
		p <- ggplot(dataForPlot1, aes(Pathway, FoldChange))
		p1 <- p + geom_boxplot(outlier.colour = NA)
		p1<-p1 + geom_point(position = position_jitter(width = 0.2),size=I(0.7))
		p2<-p1+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=20),axis.text.y=element_text(size=20))+geom_hline(yintercept = 0)+coord_flip()
		
		keggFCFigureFileName<-paste0(keggOutFileName,"/pathways_foldchange_boxplot.pdf")
		pdf(keggFCFigureFileName,width=10,height=max(7,length(unique(dataForPlot1$Pathway))*0.3))
		print(p2)
		dev.off()
	} else {
		warning(paste0("No pathway has significant fold changes! Can't make differential pathway fold change boxplot."))
	}
	
	
#Pathway Figure
#In case of unbalance fold change. balance the color
	colPart<-col_part(c(-2.5,2.5), range(resultTableFcToGene[,1]),col=colorRampPalette(c('green','black','red'))(1024))
	png(paste0(keggOutFileName,"/","FoldChange_colorBar",".png"),width=300)
	col<-col_by_value(resultTableFcToGene,col=colPart,range=c(-3,3))
	dev.off()
	
#fold change pathways
	selectedPathways<-names(which(pathwayDiff0<=pCutPathwayDiff0))
	if (length(selectedPathways)>0) {
		for (i in selectedPathways) {
			pathway_id=row.names(keggEnrichedPathway[[1]])[which(keggEnrichedPathway[[1]][,1]==i)]
#	pathway_id="04610"
#	pathway_id="04614"
#	pathway_id="04370"
#	pathway_id="03010"
#	pathway_id="04672"
#			download_KEGGfile(pathway_id)
			result_name=paste0(keggOutFileName,"/",species,"_",pathway_id,"_fc",".png")
			temp<-plot_pathway(resultTableFcToGene,type="bg",bg_col=col,text_col="white",magnify=1.4,species=species,pathway_id=pathway_id,genes_kept="abs",result_name=result_name)
		}
	} else {
		warning(paste0("No pathway has significant fold changes! Can't make differential pathway maps."))
	}
	
#differential genes enriched pathways
	sigPathwayInd=which(keggSigGeneEnrichedPathway[[1]][,pValuePathwayCol]<=pCutPathway)
	if (length(sigPathwayInd)>0) {
		for (pathway_id in row.names(keggSigGeneEnrichedPathway[[1]])[sigPathwayInd]) {
#			pathway_id=i
#			download_KEGGfile(pathway_id)
			result_name=paste0(keggOutFileName,"/",species,"_",pathway_id,"_fc",".png")
			temp<-plot_pathway(resultTableFcToGene,type="bg",bg_col=col,text_col="white",magnify=1.4,species=species,pathway_id=pathway_id,genes_kept="abs",result_name=result_name)
		}
	}
	
	
	##########################################
#output pathway information table
	##########################################
	
#temp<-find_enriched_pathway(genesEnz,species=species,returned_genenumber = 5,returned_pvalue=1,returned_adjpvalue = 1)
#	keggResultOut<-data.frame(PathwayId=row.names(keggEnrichedPathway[[1]]),stringsAsFactors=FALSE)
#	keggResultOut$PathwayId<-row.names(keggEnrichedPathway[[1]])
#	keggResultOut$PathwayId<-row.names(keggResultOut)
#	keggResultOut<-keggEnrichedPathway[[1]][,c("Gene_Found","Gene_Pathway")]
	
	keggResultOut<-data.frame(Pathway_Id=row.names(keggEnrichedPathway[[1]]),stringsAsFactors=FALSE)
	keggResultOut$Pathway_Name<-keggEnrichedPathway[[1]][,c("Pathway_Name")]
	keggResultOut$Pathway_SignificantDiffGene<-keggSigGeneEnrichedPathway[[1]][keggResultOut$Pathway_Id,c("Gene_Found")]
	keggResultOut$Pathway_DataGene<-keggEnrichedPathway[[1]][,c("Gene_Found")]
	keggResultOut$Pathway_AllGene<-keggEnrichedPathway[[1]][,c("Gene_Pathway")]
	keggResultOut$Pathway_Log2FoldChange<-pathwayMedian[keggResultOut$Pathway_Name]
	keggResultOut$Pathway_Log2FoldChangePValue<-pathwayDiff0[keggResultOut$Pathway_Name]
	keggResultOut$Pathway_SignificantDiffGenePValue<-keggSigGeneEnrichedPathway[[1]][keggResultOut$Pathway_Id,c("pvalue")]
	keggResultOut$Pathway_SignificantDiffGeneAdjPValue<-keggSigGeneEnrichedPathway[[1]][keggResultOut$Pathway_Id,c("pvalueAdj")]
	
	keggTableFileName<-paste0(keggOutFileName,".csv")
	write.csv(keggResultOut,keggTableFileName,row.names=FALSE)
	
}









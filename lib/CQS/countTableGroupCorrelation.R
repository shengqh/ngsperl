library(heatmap3)
library(DESeq2)  
library(RColorBrewer)
library(colorRamps)

countTableFileList<-parSampleFile1
groupFileList<-parSampleFile2
fixColorRange<-TRUE

if(exists("useGreenRedColorInHCA") && useGreenRedColorInHCA){
  hmcols <- colorRampPalette(c("green", "black", "red"))(256)
}else{
  hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
}

if(exists("usePearsonInHCA") && usePearsonInHCA){
  distf <- function(x) as.dist(1 - cor(t(x), use = "pa"))
}else{
  distf <- dist
}

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")

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
drawPCA<-function(filename, rldmatrix, showLabelInPCA, conditionColors){
  genecount<-nrow(rldmatrix)
  if(genecount > 2){
    cat("saving PCA to ", filename, "\n")
    png(filename=filename, width=3000, height=3000, res=300) # 10 X 10 inches
    #pdf(filename, width=10, height=10)
    pca<-prcomp(t(rldmatrix))
    supca<-summary(pca)$importance
    pcadata<-data.frame(pca$x)
    pcalabs=paste0(colnames(pcadata), "(", round(supca[2,] * 100), "%)")
    pcadata["sample"]<-row.names(pcadata)
    
    if(showLabelInPCA){
      g <- ggplot(pcadata, aes(x=PC1, y=PC2, label=sample)) + 
        geom_text(vjust=-0.6, size=4) +
        geom_point(col=conditionColors, size=4) + 
        scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) +
        scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) + 
        geom_hline(aes(yintercept=0), size=.2) + 
        geom_vline(aes(xintercept=0), size=.2) + 
        xlab(pcalabs[1]) + ylab(pcalabs[2])
    }else{
      g <- ggplot(pcadata, aes(x=PC1, y=PC2)) + 
        geom_point(col=conditionColors, size=4) + 
        labs(color = "Group") +
        scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) + 
        scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) + 
        geom_hline(aes(yintercept=0), size=.2) + 
        geom_vline(aes(xintercept=0), size=.2) +
        xlab(pcalabs[1]) + ylab(pcalabs[2]) + 
        theme(legend.position="top")
    }
    
    print(g)
    dev.off()
  }
}

##Solving node stack overflow problem end###

countTableFileAll<-read.delim(countTableFileList,header=F,as.is=T)
for (i in 1:nrow(countTableFileAll)) {
	countTableFile<-countTableFileAll[i,1]
	print(paste0("Reading ",countTableFile))
	
	if (grepl(".csv$",countTableFile)) {
		count<-read.csv(countTableFile,header=T,row.names=1,as.is=T,check.names=FALSE)
	} else {
		count<-read.delim(countTableFile,header=T,row.names=1,as.is=T,check.names=FALSE)
	}
	if (nrow(count)==0) {
		next;
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
	
	dds=DESeqDataSetFromMatrix(countData = countNum, colData = as.data.frame(rep(1,ncol(countNum))),design = ~1)

	vsdres<-try(temp<-DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
	if (class(vsdres) == "try-error") {
		message=paste0("Warning: varianceStabilizingTransformation function can't run.\n",as.character(vsdres))
		warning(message)
		writeLines(message,paste0(countTableFile,".vsd.warning"))
		next;
	}
	countNumVsd<-assay(temp)
	colnames(countNumVsd)<-colnames(countNum)
	
	#heatmap
	margin=c(max(9,max(nchar(colnames(countNumVsd)))/2), 5)
	#margin=c(min(10,max(nchar(colnames(countNumVsd)))/2),min(10,max(nchar(row.names(countNumVsd)))/2))
	
	countHT<-countNumVsd
	if(exists("top25cvInHCA") && top25cvInHCA){
    CV <- function(x){
      (sd(x)/mean(x))*100
    }
    cvs <- apply(countNumVsd,1,CV)
    countHT<-countNumVsd[cvs>=quantile(cvs)[4],]
	}
	
  if (groupFileList!="") {
    sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)
    colors=primary.colors(length(unique(sampleToGroup$V2)))
    conditionColors<-as.matrix(data.frame(Group=primary.colors(length(unique(sampleToGroup$V2)))[as.factor(sampleToGroup$V2)]))
  }else{
    conditionColors=NA
  }
  
  width=max(2000, 50 * ncol(countHT))
  print("Drawing heatmap for all samples.")
  png(paste0(countTableFile,".heatmap.png"),width=width,height=width,res=300)
  
  if(nrow(countHT) < 20){
    heatmap3(countHT,distfun=distf,margin=margin,balanceColor=TRUE,useRaster=FALSE,col=hmcols, ColSideColors=conditionColors)
  }else{
    heatmap3(countHT,distfun=distf,margin=margin,balanceColor=TRUE,useRaster=FALSE,showRowDendro=F,labRow="",col=hmcols, ColSideColors=conditionColors)
  }
  dev.off()
	
  print("Drawing PCA for all samples.")
  drawPCA(paste0(countTableFile,".PCA.png"), countHT, 1, conditionColors)
	
  print("Doing correlation analysis ...")
	#correlation distribution
	countNumCor<-cor(countNumVsd,use="pa",method="sp")
	
	colAll<-colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
	if (min(countNumCor,na.rm=T)<0) {
		colAllLabel<-c(-1,0,1)
		if (fixColorRange) {
			col<-col_part(data_all=c(-1,1),data_part=countNumCor,col=colAll)
		} else {
			col<-colAll
		}
	} else {
		colAllLabel<-c(0,0.5,1)
		if (fixColorRange) {
			col<-col_part(data_all=c(0,1),data_part=countNumCor,col=colAll)
		} else {
			col<-colAll
		}
	}
	
	legendfun<-function(x) {
		par(mar = c(5, 1, 1, 1));
		image(x=1:length(colAll),y=1,z=matrix(1:length(colAll),ncol=1),xlab="",xaxt="n",yaxt="n",col=colAll);
		axis(1,at=c(1,length(colAll)/2,length(colAll)),labels=colAllLabel)
	}
	
	png(paste0(countTableFile,".Correlation.png"),width=width,height=width,res=300)
	labRow=NULL
	margin=c(min(10,max(nchar(colnames(countNumCor)))/2),min(10,max(nchar(row.names(countNumCor)))/2))
	heatmap3(countNumCor[nrow(countNumCor):1,],scale="none",balanceColor=T,labRow=labRow,margin=margin,Rowv=NA,Colv=NA,col=col,legendfun=legendfun)
	dev.off()
	if (ncol(countNumCor)>3) {
		png(paste0(countTableFile,".Correlation.Cluster.png"),width=width,height=width,res=300)
		heatmap3(countNumCor,scale="none",balanceColor=T,labRow=labRow,margin=margin,col=col,legendfun=legendfun)
		dev.off()
	}
	
	if (groupFileList!="") {
    cexColGroup<-1
		sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)
		#keep the groups with samples in the count table
		sampleToGroup<-sampleToGroup[which(sampleToGroup[,1] %in% colnames(countNumVsd)),]
		countNumVsdGroup<-mergeTableBySampleGroup(countNumVsd,sampleToGroup)
		
		conditionColors=
		drawPCA(countTableFile, countHT, 1, designData, conditionColors)
		
		
		#heatmap
		margin=c(min(10,max(nchar(colnames(countNumVsdGroup)))/1.5),min(10,max(nchar(row.names(countNumVsdGroup)))/2))
		png(paste0(countTableFile,".Group.heatmap.png"),width=2000,height=2000,res=300)
		if(nrow(countNumVsdGroup) < 20){
			heatmap3(countNumVsdGroup,distfun=dist,margin=margin,balanceColor=TRUE,useRaster=FALSE,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),cexCol=cexColGroup)
		}else{
			margin[2]<-5
			heatmap3(countNumVsdGroup,distfun=dist,margin=margin,balanceColor=TRUE,useRaster=FALSE,labRow="",col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),cexCol=cexColGroup)
		}
		dev.off()
		
		#correlation distribution
		countNumCor<-cor(countNumVsdGroup,use="pa",method="sp")
		margin=c(min(10,max(nchar(colnames(countNumCor)))/1.5),min(10,max(nchar(row.names(countNumCor)))/1.5))
		
		colAll<-colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
#		if (min(countNumCor,na.rm=T)<0) {
#			colAllLabel<-c(-1,0,1)
#			if (fixColorRange) {
#				col<-col_part(data_all=c(-1,1),data_part=countNumCor,col=colAll)
#			} else {
#				col<-colAll
#			}
#		} else {
#			colAllLabel<-c(0,0.5,1)
#			if (fixColorRange) {
#				col<-col_part(data_all=c(0,1),data_part=countNumCor,col=colAll)
#			} else {
#				col<-colAll
#			}
#		}
		colAllLabel<-c(0,0.5,1)
		countNumCor[countNumCor<0]<-0
		if (fixColorRange) {
			col<-col_part(data_all=c(0,1),data_part=countNumCor,col=colAll)
		} else {
			col<-colAll
		}
		
		legendfun<-function(x) {
			par(mar = c(5, 1, 1, 1));
			image(x=1:length(colAll),y=1,z=matrix(1:length(colAll),ncol=1),xlab="",xaxt="n",yaxt="n",col=colAll);
			axis(1,at=c(1,length(colAll)/2,length(colAll)),labels=colAllLabel)
		}
		
		
		png(paste0(countTableFile,".Group.Correlation.png"),width=2000,height=2000,res=300)
		heatmap3(countNumCor[nrow(countNumCor):1,],scale="none",balanceColor=T,margin=margin,Rowv=NA,Colv=NA,col=col,legendfun=legendfun,cexCol=cexColGroup,cexRow=cexColGroup)
		dev.off()
		if (ncol(countNumCor)<=3 | any(is.na(cor(countNumCor,use="pa")))) {
			saveInError(paste0("Less than 3 samples. Can't do correlation analysis for group table for ",countTableFile),fileSuffix = paste0(Sys.Date(),".warning"))
		} else {
			png(paste0(countTableFile,".Group.Correlation.Cluster.png"),width=2000,height=2000,res=300)
			heatmap3(countNumCor,scale="none",balanceColor=T,margin=margin,col=col,legendfun=legendfun,cexCol=cexColGroup,cexRow=cexColGroup)
			dev.off()
		}
	}
}




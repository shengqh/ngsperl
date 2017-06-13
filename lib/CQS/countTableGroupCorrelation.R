options(bitmapType='cairo')

library(heatmap3)
library(DESeq2)  
library(RColorBrewer)
library(colorRamps)

countTableFileList<-parSampleFile1
groupFileList<-parSampleFile2
comparisonFileList<-parSampleFile3
fixColorRange<-TRUE

totalCountFile<-parFile3


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

if(!exists("showLabelInPCA")){
  showLabelInPCA<-TRUE
}

if(!exists("suffix")){
  suffix<-""
}

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
drawPCA<-function(filename, rldmatrix, showLabelInPCA, groups, groupColors){
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
        geom_hline(aes(yintercept=0), size=.2) + 
        geom_vline(aes(xintercept=0), size=.2) + 
        xlab(pcalabs[1]) + ylab(pcalabs[2])
    
    print(g)
    dev.off()
  }
}

##Solving node stack overflow problem end###

#correlation without 0
corWithout<-function (x, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman"),ignorValue=0) 
{
	na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
					"everything", "na.or.complete"))
	if (is.na(na.method)) 
		stop("invalid 'use' argument")
	method <- match.arg(method)
	if (is.data.frame(y)) 
		y <- as.matrix(y)
	if (is.data.frame(x)) 
		x <- as.matrix(x)
	if (!is.matrix(x) && is.null(y)) 
		stop("supply both 'x' and 'y' or a matrix-like 'x'")
	if (!(is.numeric(x) || is.logical(x))) 
		stop("'x' must be numeric")
	stopifnot(is.atomic(x))
	if (!is.null(y)) {
		if (!(is.numeric(y) || is.logical(y))) 
			stop("'y' must be numeric")
		stopifnot(is.atomic(y))
	}
	Rank <- function(u) {
		if (length(u) == 0L) 
			u
		else if (is.matrix(u)) {
			if (nrow(u) > 1L) 
				apply(u, 2L, rank, na.last = "keep")
			else row(u)
		}
		else rank(u, na.last = "keep")
	}
	
	#codes for ignor value
	if (!is.null(ignorValue)) {
		x[x==ignorValue]<-NA
		if (!is.null(y)) {
			y[y==ignorValue]<-NA
		}
	}
	
	if (method == "pearson") 
		.Call(stats:::C_cor, x, y, na.method, FALSE)
	else if (na.method %in% c(2L, 5L)) {
		if (is.null(y)) {
			.Call(stats:::C_cor, Rank(na.omit(x)), NULL, na.method, method == 
							"kendall")
		}
		else {
			nas <- attr(na.omit(cbind(x, y)), "na.action")
			dropNA <- function(x, nas) {
				if (length(nas)) {
					if (is.matrix(x)) 
						x[-nas, , drop = FALSE]
					else x[-nas]
				}
				else x
			}
			.Call(stats:::C_cor, Rank(dropNA(x, nas)), Rank(dropNA(y, 
									nas)), na.method, method == "kendall")
		}
	}
	else if (na.method != 3L) {
		x <- Rank(x)
		if (!is.null(y)) 
			y <- Rank(y)
		.Call(stats:::C_cor, x, y, na.method, method == "kendall")
	}
	else {
		if (is.null(y)) {
			ncy <- ncx <- ncol(x)
			if (ncx == 0) 
				stop("'x' is empty")
			r <- matrix(0, nrow = ncx, ncol = ncy)
			for (i in seq_len(ncx)) {
				for (j in seq_len(i)) {
					x2 <- x[, i]
					y2 <- x[, j]
					ok <- complete.cases(x2, y2)
					x2 <- rank(x2[ok])
					y2 <- rank(y2[ok])
					r[i, j] <- if (any(ok)) 
								.Call(stats:::C_cor, x2, y2, 1L, method == "kendall")
							else NA
				}
			}
			r <- r + t(r) - diag(diag(r))
			rownames(r) <- colnames(x)
			colnames(r) <- colnames(x)
			r
		}
		else {
			if (length(x) == 0L || length(y) == 0L) 
				stop("both 'x' and 'y' must be non-empty")
			matrix_result <- is.matrix(x) || is.matrix(y)
			if (!is.matrix(x)) 
				x <- matrix(x, ncol = 1L)
			if (!is.matrix(y)) 
				y <- matrix(y, ncol = 1L)
			ncx <- ncol(x)
			ncy <- ncol(y)
			r <- matrix(0, nrow = ncx, ncol = ncy)
			for (i in seq_len(ncx)) {
				for (j in seq_len(ncy)) {
					x2 <- x[, i]
					y2 <- y[, j]
					ok <- complete.cases(x2, y2)
					x2 <- rank(x2[ok])
					y2 <- rank(y2[ok])
					r[i, j] <- if (any(ok)) 
								.Call(stats:::C_cor, x2, y2, 1L, method == "kendall")
							else NA
				}
			}
			rownames(r) <- colnames(x)
			colnames(r) <- colnames(y)
			if (matrix_result) 
				r
			else drop(r)
		}
	}
}


#start work:
countTableFileAll<-read.delim(countTableFileList,header=F,as.is=T)
i<-1
for (i in 1:nrow(countTableFileAll)) {
  countTableFile<-countTableFileAll[i,1]
  
  if(outputDirectory==""){
    outputFilePrefix=countTableFile
  }else{
    bname=basename(countTableFile)
    dpath=dirname(countTableFile)
    dname=basename(dpath)
    while(dname=="result"){
      dpath=dirname(dpath)
      dname=basename(dpath)
    }
    outputFilePrefix=paste0(outputDirectory, "/", dname, "." , bname)
  }
  
  countTableTitle<-countTableFileAll[i,2]
  
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
  #remove samples with total reads 0
  countNum<-countNum[,which(colSums(countNum,na.rm=T)>0)]
  
  if (groupFileList!="") {
    sampleToGroup<-getSampleInGroup(groupFileList, colnames(countNum), comparisonFileList, countTableTitle, useLeastGroups)
    write.table(sampleToGroup, paste0(outputFilePrefix,suffix,".correlation.groups"),col.names=F, row.names=F, quote=F, sep="\t")
  }
  
  #filter reads/genes by parameter
  countNum<-filterCountTable(countNum,groupFileList=groupFileList,minMedian=minMedian,minMedianInGroup=minMedianInGroup)
  
  #normlize by total reads or VSD
  if (totalCountFile!="") { #normlize with total count *10^6
	  totalCount<-read.csv(totalCountFile,header=T,as.is=T,row.names=1,check.names=FALSE)
	  totalCount<-unlist(totalCount[totalCountKey,])
	  countNumVsd<-10^6*t(t(countNum)/totalCount[colnames(countNum)])
	  write.table(countNumVsd, paste0(outputFilePrefix,suffix,".RPM.txt"),col.names=NA, quote=F, sep="\t")
	  
	  ylab<-"Mapped Reads per Million"
  } else {
	  dds=DESeqDataSetFromMatrix(countData = countNum, colData = as.data.frame(rep(1,ncol(countNum))),design = ~1)
	  dds<-try(myEstimateSizeFactors(dds))
	  vsdres<-try(temp<-DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
	  if (class(vsdres) == "try-error") {
		  message=paste0("Warning: varianceStabilizingTransformation function failed.\n",as.character(vsdres))
		  warning(message)
		  writeLines(message,paste0(outputFilePrefix,suffix,".vsd.warning"))
		  next;
	  }
	  countNumVsd<-assay(temp)
	  colnames(countNumVsd)<-colnames(countNum)
	  write.table(countNumVsd, paste0(outputFilePrefix,suffix,".vsd.txt"),col.names=NA, quote=F, sep="\t")
	  
	  ylab<-"VSD"
  }
  
#  #To be removed
#  dds=DESeqDataSetFromMatrix(countData = countNum, colData = as.data.frame(rep(1,ncol(countNum))),design = ~1)
#  dds<-try(myEstimateSizeFactors(dds))
#  vsdres<-try(temp<-DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
#  if (class(vsdres) == "try-error") {
#    message=paste0("Warning: varianceStabilizingTransformation function failed.\n",as.character(vsdres))
#    warning(message)
#    writeLines(message,paste0(outputFilePrefix,suffix,".vsd.warning"))
#    next;
#  }
#  countNumVsd<-assay(temp)
#  colnames(countNumVsd)<-colnames(countNum)
#  write.table(countNumVsd, paste0(outputFilePrefix,suffix,".vsd.txt"),col.names=NA, quote=F, sep="\t")
#  #end to be removed
  
  
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
    groups<-sampleToGroup$V2
    colors<-primary.colors(length(unique(groups)))
    conditionColors<-as.matrix(data.frame(Group=colors[groups]))
  }else{
    groups<-NA
    colors<-NA
    conditionColors<-NA
  }
  
  print("Drawing PCA for all samples.")
  drawPCA(paste0(outputFilePrefix,suffix,".PCA.png"), countHT, showLabelInPCA, groups, colors)
  
  width=max(2000, 50 * ncol(countHT))
  print("Drawing heatmap for all samples.")
  png(paste0(outputFilePrefix,suffix,".heatmap.png"),width=width,height=width,res=300)
  
  if(nrow(countHT) < 20){
    if(!is.na(conditionColors[1])){
      heatmap3(countHT,distfun=distf,margin=margin,balanceColor=TRUE,useRaster=FALSE,col=hmcols, ColSideColors=conditionColors)
    }else{
      heatmap3(countHT,distfun=distf,margin=margin,balanceColor=TRUE,useRaster=FALSE,col=hmcols)
    }
  }else{
    if(!is.na(conditionColors[1])){
      heatmap3(countHT,distfun=distf,margin=margin,balanceColor=TRUE,useRaster=FALSE,showRowDendro=F,labRow="",col=hmcols, ColSideColors=conditionColors)
    }else{
      heatmap3(countHT,distfun=distf,margin=margin,balanceColor=TRUE,useRaster=FALSE,showRowDendro=F,labRow="",col=hmcols)
    }
  }
  dev.off()
  
  print("Doing correlation analysis ...")
  #correlation distribution
  countNumCor<-corWithout(countNumVsd,use="pa",method="sp")
  
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
  if (groupFileList!="") {
	  legendfun<-function() showLegend(legend=unique(groups),col=unique(conditionColors[,1]))
  }
  
  png(paste0(outputFilePrefix,suffix,".Correlation.png"),width=width,height=width,res=300)
  labRow=NULL
  margin=c(min(10,max(nchar(colnames(countNumCor)))/2),min(10,max(nchar(row.names(countNumCor)))/2))
  if(all(!is.na(conditionColors))) { #has group information
	  heatmap3(countNumCor[nrow(countNumCor):1,],scale="none",balanceColor=T,labRow=labRow,margin=margin,Rowv=NA,Colv=NA,col=col,legendfun=legendfun,ColSideColors=conditionColors)
  }else{
	  heatmap3(countNumCor[nrow(countNumCor):1,],scale="none",balanceColor=T,labRow=labRow,margin=margin,Rowv=NA,Colv=NA,col=col,legendfun=legendfun)
  }
  dev.off()
  if (ncol(countNumCor)>3) {
    png(paste0(outputFilePrefix,suffix,".Correlation.Cluster.png"),width=width,height=width,res=300)
	if(!is.na(conditionColors[1,1])){
		heatmap3(countNumCor,scale="none",balanceColor=T,labRow=labRow,margin=margin,col=col,legendfun=legendfun,ColSideColors=conditionColors)
	}else{
		heatmap3(countNumCor,scale="none",balanceColor=T,labRow=labRow,margin=margin,col=col,legendfun=legendfun)
	}
    dev.off()
  }
  
  if (groupFileList!="") {
    cexColGroup<-1
    if(length(unique(sampleToGroup$V2)) < 3){
      saveInError(paste0("Less than 3 groups. Can't do correlation analysis for group table for ",countTableFile),fileSuffix = paste0(suffix,Sys.Date(),".warning"))
      next
    }

    countNumVsdGroup<-mergeTableBySampleGroup(countNumVsd,sampleToGroup)
    
    #heatmap
    margin=c(min(10,max(nchar(colnames(countNumVsdGroup)))/1.5),min(10,max(nchar(row.names(countNumVsdGroup)))/2))
    png(paste0(outputFilePrefix,suffix,".Group.heatmap.png"),width=2000,height=2000,res=300)
    if(nrow(countNumVsdGroup) < 20){
      heatmap3(countNumVsdGroup,distfun=dist,margin=margin,balanceColor=TRUE,useRaster=FALSE,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),cexCol=cexColGroup)
    }else{
      margin[2]<-5
      heatmap3(countNumVsdGroup,distfun=dist,margin=margin,balanceColor=TRUE,useRaster=FALSE,labRow="",col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),cexCol=cexColGroup)
    }
    dev.off()
    
    #correlation distribution
    countNumCor<-corWithout(countNumVsdGroup,use="pa",method="sp")
    margin=c(min(10,max(nchar(colnames(countNumCor)))/1.5),min(10,max(nchar(row.names(countNumCor)))/1.5))
    
    colAll<-colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
#   if (min(countNumCor,na.rm=T)<0) {
#     colAllLabel<-c(-1,0,1)
#     if (fixColorRange) {
#       col<-col_part(data_all=c(-1,1),data_part=countNumCor,col=colAll)
#     } else {
#       col<-colAll
#     }
#   } else {
#     colAllLabel<-c(0,0.5,1)
#     if (fixColorRange) {
#       col<-col_part(data_all=c(0,1),data_part=countNumCor,col=colAll)
#     } else {
#       col<-colAll
#     }
#   }
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
    
    png(paste0(outputFilePrefix,suffix,".Group.Correlation.png"),width=2000,height=2000,res=300)
    heatmap3(countNumCor[nrow(countNumCor):1,],scale="none",balanceColor=T,margin=margin,Rowv=NA,Colv=NA,col=col,legendfun=legendfun,cexCol=cexColGroup,cexRow=cexColGroup)
    dev.off()
    if (ncol(countNumCor)<=3 | any(is.na(corWithout(countNumCor,use="pa")))) {
      saveInError(paste0("Less than 3 samples. Can't do correlation analysis for group table for ",countTableFile),fileSuffix = paste0(outputFilePrefix,suffix,Sys.Date(),".warning"))
    } else {
      png(paste0(outputFilePrefix,suffix,".Group.Correlation.Cluster.png"),width=2000,height=2000,res=300)
      heatmap3(countNumCor,scale="none",balanceColor=T,margin=margin,col=col,legendfun=legendfun,cexCol=cexColGroup,cexRow=cexColGroup)
      dev.off()
    }
  }
}

countTableFileList<-parSampleFile1
groupFileList<-parSampleFile2
fixColorRange<-TRUE

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")

library(heatmap3)
library(DESeq2)  
library(RColorBrewer)

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
	
	colClass<-sapply(count, class)
	countNotNumIndex<-which(colnames(count) == "Feature_length" | (colClass!="numeric" & colClass!="integer"))
	if (length(countNotNumIndex)==0) {
		countNotNumIndex<-0;
	} else {
		countNotNumIndex<-max(countNotNumIndex)
	}
	countNum<-count[,c((countNotNumIndex+1):ncol(count))]
	countNum<-round(countNum,0)
	#remove genes with total reads 0
	countNum<-countNum[which(rowSums(countNum,na.rm=T)>0),]
	
	dds=DESeqDataSetFromMatrix(countData = countNum, colData = as.data.frame(rep(1,ncol(countNum))),design = ~1)
	temp<-DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
	countNumVsd<-assay(temp)
	colnames(countNumVsd)<-colnames(countNum)
	
	#heatmap
	margin=c(max(9,max(nchar(colnames(countNumVsd)))/2), 5)
	#margin=c(min(10,max(nchar(colnames(countNumVsd)))/2),min(10,max(nchar(row.names(countNumVsd)))/2))
	png(paste0(countTableFile,".heatmap.png"),width=2000,height=2000,res=300)
	if(nrow(countNum) < 20){
	  heatmap3(countNumVsd,dist=dist,margin=margin,balanceColor=TRUE,useRaster=FALSE,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
	}else{
	  heatmap3(countNumVsd,dist=dist,margin=margin,balanceColor=TRUE,useRaster=FALSE,labRow="",col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
	}
	dev.off()
	
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
	
	png(paste0(countTableFile,".Correlation.png"),width=2000,height=2000,res=300)
	if(nrow(countNum) < 20){
	  margin=c(min(10,max(nchar(colnames(countNumCor)))/2),min(10,max(nchar(row.names(countNumCor)))/2))
	  heatmap3(countNumCor[nrow(countNumCor):1,],scale="none",balanceColor=T,margin=margin,Rowv=NA,Colv=NA,col=col,legendfun=legendfun)
	}else{
	  margin=c(min(10,max(nchar(colnames(countNumCor)))/2),5)
	  heatmap3(countNumCor[nrow(countNumCor):1,],scale="none",balanceColor=T,labRow="",margin=margin,Rowv=NA,Colv=NA,col=col,legendfun=legendfun)
	}
	dev.off()
	if (ncol(countNumCor)>3) {
		png(paste0(countTableFile,".Correlation.Cluster.png"),width=2000,height=2000,res=300)
		heatmap3(countNumCor,scale="none",balanceColor=T,margin=margin,col=col,legendfun=legendfun)
		dev.off()
	}
	
	if (groupFileList!="") {
		sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)
		#keep the groups with samples in the count table
		sampleToGroup<-sampleToGroup[which(sampleToGroup[,1] %in% colnames(countNumVsd)),]
		countNumVsdGroup<-mergeTableBySampleGroup(countNumVsd,sampleToGroup)
		
		#heatmap
		margin=c(min(10,max(nchar(colnames(countNumVsdGroup)))/2),min(10,max(nchar(row.names(countNumVsdGroup)))/2))
		png(paste0(countTableFile,".Group.heatmap.png"),width=2000,height=2000,res=300)
		heatmap3(countNumVsdGroup,dist=dist,margin=margin,balanceColor=TRUE,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
		dev.off()
		
		#correlation distribution
		countNumCor<-cor(countNumVsdGroup,use="pa",method="sp")
		margin=c(min(10,max(nchar(colnames(countNumCor)))/2),min(10,max(nchar(row.names(countNumCor)))/2))
		
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
		
		
		png(paste0(countTableFile,".Group.Correlation.png"),width=2000,height=2000,res=300)
		heatmap3(countNumCor[nrow(countNumCor):1,],scale="none",balanceColor=T,margin=margin,Rowv=NA,Colv=NA,col=col,legendfun=legendfun)
		dev.off()
		if (ncol(countNumCor)<=3 | any(is.na(cor(countNumCor,use="pa")))) {
			saveInError(paste0("Can't do correlation analysis for group table for ",countTableFile))
		} else {
			png(paste0(countTableFile,".Group.Correlation.Cluster.png"),width=2000,height=2000,res=300)
			heatmap3(countNumCor,scale="none",balanceColor=T,margin=margin,col=col,legendfun=legendfun)
			dev.off()
		}
	}
}




countTableFileList<-parSampleFile1
groupFileList<-parSampleFile2

library(heatmap3)
library(DESeq2)  
library(RColorBrewer)
library(ggplot2)


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

text2Color<-function(x) {
  colNum=max(length(unique(x)),3)
  col=RColorBrewer::brewer.pal(colNum,"Set1")[1:length(unique(x))]
  names(col)<-unique(x)
  colOut<-col[x]
  return(list(color=colOut,legend=col))
}
drawPCA<-function(prefix, rldmatrix, showLabelInPCA, conditionColors){
  filename<-paste0(prefix, ".PCA.png")
  genecount<-nrow(rldmatrix)
  if(genecount > 2){
    cat("saving PCA to ", filename, "\n")
    png(filename=filename, width=2000, height=2000, res=300)
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

panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  #	r <- abs(cor(x, y))
  r <- (cor(x, y,use="pa",method="spearman"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y,method="spearman")
  Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
  text(0.5, 0.25, paste("r=",txt))
  text(.5, .75, Signif)
}

panel.smooth<-function (x, y, col = "blue", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}


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
  
  #Group
  sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)
  #keep the groups with samples in the count table
  sampleToGroup<-sampleToGroup[which(sampleToGroup[,1] %in% colnames(countNumVsd)),]
  countNumVsdOrdered<-countNumVsd[,sampleToGroup[,1]]
  groupColor<-text2Color(sampleToGroup[,2])$color
  
  #heatmap
  margins=c(max(9,max(nchar(colnames(countNumVsd)))/2), 5)
  if (nrow(countNumVsdOrdered)>10) {
	  labRow=""
  } else {
	  labRow=row.names(countNumVsdOrdered)
  }
  png(paste0(countTableFile,".heatmap.png"),width=2000,height=2000,res=300)
  heatmap3(countNumVsdOrdered,ColSideColors = groupColor,ColSideLabs="Group",labRow=labRow, dist=dist,balanceColor=TRUE,useRaster=FALSE,margins=margins)
  dev.off()
  
  #PCA
  drawPCA(countTableFile, countNumVsdOrdered, showLabelInPCA=TRUE, groupColor)
  
  #Pairs correlation
  swidth=max(2000, ncol(countNumVsdOrdered) * 400)
  png(paste0(countTableFile,".pairsCorrelation.png"),width=swidth,height=swidth,res=300)
  pairs(countNumVsd,lower.panel=panel.smooth, upper.panel=panel.cor)
  dev.off()
}




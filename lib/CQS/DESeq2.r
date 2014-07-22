##predefined_condition_begin
#  setwd("H:/shengquanhu/projects/chenxi/20131017_chenxi_rnaseq_smad4/deseq2/result")  
#  data<-read.table("H:/shengquanhu/projects/chenxi/20131017_chenxi_rnaseq_smad4/genetable/result/smad4_gene.count",row.names=1, header=T, check.names=F)
#  pairs=list(
#  	"KO_vs_WT" = list("WT" = c("2288-RDB-81","2288-RDB-83","2288-RDB-85"), 
#                     "KO" = c("2288-RDB-82","2288-RDB-84","2288-RDB-86"), 
#                     "paired" = c("S1", "S2", "S3")), 
#    "KO_vs_WT_unpair" = list("WT" = c("2288-RDB-81","2288-RDB-83","2288-RDB-85"), "KO" = c("2288-RDB-82","2288-RDB-84","2288-RDB-86")) 
# )
##predefined_condition_end

library("DESeq2")
library("heatmap3")
library("lattice")
library("reshape")
library("ggplot2")
library("grid")

hmcols <- colorRampPalette(c("green", "black", "red"))(256)

countData<-data
index<-1
indecies<-c()
while(! is.numeric(countData[1,1])){
  countData<-countData[,c(2:ncol(countData))]
  indecies<-c(indecies, index)
  index<-index+1
}
countData[is.na(countData)] <- 0
countData<-round(countData)

pairnames=names(pairs)
pairname=pairnames[1]

pairedspearman<-list()

for(pairname in pairnames){
  str(pairname)
  gs=pairs[[pairname]]
  gnames=names(gs)
  if(length(gnames) > 2){
    ispaired<-TRUE
    pairsamplenames<-unlist(gs[3])
    cat("Paired data!\n")
  }else{
    ispaired<-FALSE
    cat("Not paired data!\n")
  }
  
  gnames<-gnames[1:2]
  g1name=gnames[1]
  g2name=gnames[2]
  g1=gs[[g1name]]
  g2=gs[[g2name]]
  c1=countData[,colnames(countData) %in% g1,drop=F]
  c2=countData[,colnames(countData) %in% g2,drop=F]
  
  if(ncol(c1) != length(g1)){
    warning(paste0("There are only ", ncol(c1), " samples in group ", g1name, " but ", length(g1), " required!"))
    next
  }
  
  if(ncol(c2) != length(g2)){
    warning(paste0("There are only ", ncol(c2), " samples in group ", g2name, " but ", length(g2), " required!"))
    next
  }
  c1<-c1[,g1]
  c2<-c2[,g2]
  
  if(ispaired){
    if(ncol(c2) != ncol(c1)){
      warning(paste0("Data not paired, there are ", ncol(c2), " samples in group ", g2name, " but ", ncol(c1), " samples in group ", g1name))
      next
    }

    if(ncol(c2) != length(pairsamplenames)){
      warning(paste0("Name not paired with sample, there are ", ncol(c2), " samples in group ", g2name, " but ", length(pairsamplenames), " names defined"))
      next
    }
    
    spcorr<-unlist(lapply(c(1:length(g1)), function(x){
              cor(c1[,x], c2[,x],method="spearman")
            }))
            

    sptable<-data.frame(Name=pairsamplenames, Spcorr=spcorr)
    write.csv(sptable, file=paste0(pairname, "_Spearman.csv"), row.names=FALSE)
    
    lapply(c(1:length(g1)), function(x){
      log2c1<-log2(c1[,x]+1)
      log2c2<-log2(c2[,x]+1)
      png(paste0(pairname, "_Spearman_", pairsamplenames[x], ".png"), width=2000, height=2000, res=300)
      plot(log2c1, log2c2, xlab=paste0(g1[x], " [log2(Count + 1)]"), ylab=paste0(g2[x], " [log2(Count + 1)]"))
      text(3,15,paste0("SpearmanCorr=", sprintf("%0.3f", cor(c1[,x], c2[,x],method="spearman")) ))
      dev.off()
    })
    
    pairedspearman[[pairname]]<-spcorr
  }
  
  pairCountData<-cbind(c1, c2)
  notEmptyData<-apply(pairCountData, 1, max) > 0
  pairCountData<-pairCountData[notEmptyData,]
  
  if(ispaired){
    pairColData=data.frame(condition=factor(c(rep(g1name, ncol(c1)), rep(g2name, ncol(c2))), levels=gnames[1:2]), paired=factor(c(rep(pairsamplenames,2))))
    colnames(pairCountData)<-unlist(lapply(c(1:ncol(pairCountData)), function(i){paste0(pairColData$paired[i], "_", colnames(pairCountData)[i])}))
  }else{
    pairColData=data.frame(condition=factor(c(rep(g1name, ncol(c1)), rep(g2name, ncol(c2))), levels=gnames[1:2]))
  }
  rownames(pairColData)<-colnames(pairCountData)
  pairColors<-as.matrix(data.frame(Group=c(rep("red", ncol(c1)), rep("blue", ncol(c2)))))
  
  #different expression analysis
  if(ispaired){
    dds=DESeqDataSetFromMatrix(countData = pairCountData,
        colData = pairColData,
        design = ~ paired + condition)
  }else{
    dds=DESeqDataSetFromMatrix(countData = pairCountData,
        colData = pairColData,
        design = ~ condition)
  }
  
  dds <- DESeq(dds)
  res<-results(dds,cooksCutoff=FALSE)
  
  cat("DESeq2 finished.\n")
  
  select<- (!is.na(res$padj)) & (res$padj<0.05) & ((res$log2FoldChange >= 1) | (res$log2FoldChange <= -1))
  
  if(length(indecies) > 0){
    inddata<-data[notEmptyData,indecies,drop=F]
    cat(paste0(nrow(inddata)), " : ", nrow(pairCountData), " : ", nrow(res), "\n")
    tbb<-cbind(inddata, pairCountData, res)
  }else{
    tbb<-cbind(pairCountData, res)
  }
  tbbselect<-tbb[select,,drop=F]
  
  tbb<-tbb[order(tbb$padj),,drop=F]
  write.csv(as.data.frame(tbb),paste0(pairname, "_DESeq2.csv"))
  
  tbbselect<-tbbselect[order(tbbselect$padj),,drop=F]
  write.csv(as.data.frame(tbbselect),paste0(pairname, "_DESeq2_sig.csv"))
  
  #some basic graph
  dds=DESeqDataSetFromMatrix(countData = pairCountData,
      colData = pairColData,
      design = ~1)
  
  colnames(dds)<-colnames(pairCountData)
  
  #draw density graph
  rldmatrix<-as.matrix(log2(counts(dds,normalized=FALSE) + 1))
  rsdata<-melt(rldmatrix)
  colnames(rsdata)<-c("Gene", "Sample", "log2Count")
  png(filename=paste0(pairname, "_DESeq2-log2-density.png"), width=4000, height=3000, res=300)
  g<-ggplot(rsdata) + geom_density(aes(x=log2Count, colour=Sample)) + xlab("DESeq2 log2 transformed count")
  print(g)
  dev.off()
  
  #varianceStabilizingTransformation
  vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
  assayvsd<-assay(vsd)
  write.csv(assayvsd, file=paste0(pairname, "_DESeq2-vsd.csv"))
  
  vsdiqr<-apply(assayvsd, 1, IQR)
  assayvsd<-assayvsd[order(vsdiqr, decreasing=T),]
  
  rldmatrix=as.matrix(assayvsd)
  
  #draw pca graph
  png(filename=paste0(pairname, "_DESeq2-vsd-pca.png"), width=3000, height=3000, res=300)
  pca<-prcomp(t(rldmatrix))
  supca<-summary(pca)$importance
  pcadata<-data.frame(pca$x)
  pcalabs=paste0(colnames(pcadata), "(", round(supca[2,] * 100), "%)")
  g <- ggplot(pcadata, aes(x=PC1, y=PC2, label=row.names(pcadata))) + 
      geom_text(vjust=-0.6, size=4) +
      geom_point(col=pairColors, size=4) + 
      scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) +
      scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) + 
      geom_hline(aes(0), size=.2) + 
      geom_vline(aes(0), size=.2) + 
      xlab(pcalabs[1]) + ylab(pcalabs[2])
  print(g)
  dev.off()
  
  #draw heatmap
  rldselect<-rldmatrix[1:500,,drop=F]
  htfile<-paste0(pairname, "_DESeq2-vsd-heatmap.png")
  if(!file.exists(htfile)){
    if(nrow(rldselect) > 2){
      png(filename=htfile, width=3000, height =3000, res=300)
      if(ispaired){
        htColors<-rainbow(ncol(c1))
        gsColors<-as.matrix(data.frame(Group=c(rep("red", ncol(c1)), rep("blue", ncol(c2))),
                Sample=c(rep(htColors, 2))))
        heatmap3(rldselect, col = hmcols, ColSideColors = gsColors, margins=c(12,5), scale="r", dist=dist, labRow="",
            legendfun=function() showLegend(legend=paste0("Group ", gnames),col=c("red","blue"),cex=1.0,x="center"))
      }else{
        heatmap3(rldselect, col = hmcols, ColSideColors = pairColors, margins=c(12,5), scale="r", dist=dist, labRow="",
            legendfun=function() showLegend(legend=paste0("Group ", gnames),col=c("red","blue"),cex=1.0,x="center"))
      }
      dev.off()
    }
  }
}

if(length(pairedspearman) > 0){
  #draw pca graph
  png(filename=paste0("spearman.png"), width=1000 * length(pairedspearman), height=2000, res=300)
  boxplot(pairedspearman)
  dev.off()
}

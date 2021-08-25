



sciCloneInputListFile=parFile1
#sciCloneInputListFile="/gpfs23/scratch/massion_lab/BrushingUMIWES/20200730_UMI_MuTect2PON_simpleProcessing/files_CollectAllelicCounts_PrepareClonalAnalysis/result/UMIMutect2PON.sciCloneInputList.txt"

#run sciClone
library(sciClone)
minimumDepth=50
maximumClusters=10

sciCloneInputListTable=read.delim(sciCloneInputListFile,header=T,as.is=T)
sciCloneResultForClonevol=data.frame()
for (i in 1:nrow(sciCloneInputListTable)) {
  inputPatient=sciCloneInputListTable[i,1]
  outputName=paste0(inputPatient,".sciClone")
  
  inputObjOne=readRDS(sciCloneInputListTable[i,2])
  AllelicCountsTableListSciClone=inputObjOne[[1]]
  cnvTableListSciClone=inputObjOne[[2]]
  sampleNames=inputObjOne[[3]]
  
  sc = try(sciClone(vafs=AllelicCountsTableListSciClone,
                    copyNumberCalls=cnvTableListSciClone,
                    sampleNames=sampleNames,
                    minimumDepth=minimumDepth,
                    #cnCallsAreLog2=TRUE,
                    doClusteringAlongMargins=FALSE,
                    maximumClusters=maximumClusters
  ))
  if (is.null(sc)) {
    print(paste0("Not enough valid variants to run sciClone in Patient: ", inputPatient))
    next;
  }
  
  writeClusterTable(sc, paste0(outputName,".txt"))
  sc.plot1d(sc,paste0(outputName,".clusterVis1D.pdf"))
  sc.plot2d(sc,paste0(outputName,".clusterVis2D.pdf"))
  
  sciCloneResultForClonevol=rbind(sciCloneResultForClonevol,c(inputPatient,paste0(getwd(),"/",outputName,".txt")))
}

colnames(sciCloneResultForClonevol)=c("Sample","File")





#run clonevol
library(clonevol)

for (i in 1:nrow(sciCloneResultForClonevol)) {
  inputPatient=sciCloneResultForClonevol[i,1]
  cloneResultFile=sciCloneResultForClonevol[i,2]
  #cloneResultFile="6359.sciClone.txt"
  
  cloneResult=read.delim(cloneResultFile,header=T,as.is=T,check.names=F)
  
  ########################################
  #format sciClone result, export cluster column and vaf column
  ########################################
  clusterCol="cluster"
  vafCol=grep(".vaf$",colnames(cloneResult),value=TRUE)
  
  dataForVis=cloneResult[,c(clusterCol,vafCol)]
  dataForVis=dataForVis[which(!is.na(dataForVis[,clusterCol])),]
  dataForVis=dataForVis[which(dataForVis[,clusterCol]!=0),]
  dataForVis=dataForVis[order(dataForVis[,clusterCol]),]
  #remove vaf in vafCol to make it easier for viewing
  colnames(dataForVis)=gsub(".vaf","",colnames(dataForVis))
  vafCol=gsub(".vaf","",vafCol)
  
  ########################################
  #run clonevol
  ########################################
  y  = try(infer.clonal.models(
    variants = dataForVis,
    cluster.col.name =clusterCol,
    vaf.col.names = vafCol,
    cancer.initiation.model = 'monoclonal',
    subclonal.test = 'bootstrap',
    subclonal.test.model = 'non-parametric',
    num.boots = 1000,
    founding.cluster = 1,
    cluster.center = 'mean',
    ignore.clusters = NULL,
    min.cluster.vaf = 0.01,
    # min probability that CCF(clone) is non-negative
    sum.p = 0.05,
    # alpha level in confidence interval estimate for CCF(clone)
    alpha = 0.05
  ))
  
  if (class(y)=="try-error" | is.null(y)) {
    print(paste0("clonevol for ",basename(cloneResultFile)," not sucvcessful. Too few command variants between samples?"))
    next
  }
  
  y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')
  
  
  try(
  plot.clonal.models(
    y,
    # box plot parameters
    box.plot = TRUE,
    fancy.boxplot = TRUE,
#    fancy.variant.boxplot.highlight = 'is.driver',
#    fancy.variant.boxplot.highlight.shape = 21,
#    fancy.variant.boxplot.highlight.fill.color = 'red',
#    fancy.variant.boxplot.highlight.color = 'black',
#    fancy.variant.boxplot.highlight.note.col.name = 'gene',
#    fancy.variant.boxplot.highlight.note.color = 'blue',
#    fancy.variant.boxplot.highlight.note.size = 2,
    fancy.variant.boxplot.jitter.alpha = 1,
    fancy.variant.boxplot.jitter.center.color = 'grey50',
    fancy.variant.boxplot.base_size = 12,
    fancy.variant.boxplot.plot.margin = 1,
    fancy.variant.boxplot.vaf.suffix = '',
    # bell plot parameters
    clone.shape = 'bell',
    bell.event = TRUE,
    bell.event.label.color = 'blue',
    bell.event.label.angle = 60,
    clone.time.step.scale = 1,
    bell.curve.step = 2,
    # node-based consensus tree parameters
    merged.tree.plot = TRUE,
    tree.node.label.split.character = NULL,
    tree.node.shape = 'circle',
    tree.node.size = 30,
    tree.node.text.size = 0.5,
    merged.tree.node.size.scale = 1.25,
    merged.tree.node.text.size.scale = 2.5,
    merged.tree.cell.frac.ci = FALSE,
    # branch-based consensus tree parameters
    merged.tree.clone.as.branch = TRUE,
    mtcab.event.sep.char = ',',
    mtcab.branch.text.size = 1,
    mtcab.branch.width = 0.75,
    mtcab.node.size = 3,
    mtcab.node.label.size = 1,
    mtcab.node.text.size = 1.5,
    # cellular population parameters
    cell.plot = TRUE,
    num.cells = 100,
    cell.border.size = 0.25,
    cell.border.color = 'black',
    clone.grouping = 'horizontal',
    #meta-parameters
    scale.monoclonal.cell.frac = TRUE,
    show.score = FALSE,
    cell.frac.ci = TRUE,
    disable.cell.frac = FALSE,
    # output figure parameters
    #out.dir = basename(cloneResultFile),
    out.dir = ".",
    out.prefix= basename(cloneResultFile),
    out.format = 'pdf',
    overwrite.output = TRUE,
    width = 15,
    height = 5,
    # vector of width scales for each panel from left to right
    panel.widths = c(3, 4, 2, 4, 2)
  )
  )
#  break;
}






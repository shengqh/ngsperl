library(tidyverse)
require(clusterProfiler)
require(ReactomePA)
require(tidytext)
suppressPackageStartupMessages(library(ComplexHeatmap))

add_coverage<-function(dataForPlotFrame){
  dataForPlotFrame$BgCount = as.numeric(gsub("/.+","",dataForPlotFrame$BgRatio))
  dataForPlotFrame$Coverage = dataForPlotFrame$Count / dataForPlotFrame$BgCount * 100
  dataForPlotFrame
}

#' example
#' temp=c("CA1","ICAM1","IGFBP3","SERPINA5","CCL14","PAM")
#' result=enrichmentByClusterProfiler(temp,modules=c("Reactome"))
enrichmentByClusterProfiler = function(selectedProteins,
                                       universe = NULL,
                                       modules=c("Reactome","KEGG","GO.BP","GO.MF","GO.CC","WikiPathways"),
                                       organism = "hsa",
                                       fromType = "SYMBOL",
                                       toType = "ENTREZID",
                                       OrgDb = "org.Hs.eg.db",
                                       minGSSize = 3,
                                       maxGSSize = 1000,
                                       pvalueCutoff = 1,
                                       qvalueCutoff = 1,
                                       showConvertedGeneId=TRUE,
                                       showProgress = TRUE) {

  listEnrichments = list()

  if (!is.null(fromType) & !is.null(toType)) {
    eg = bitr(
      selectedProteins,
      fromType = fromType,
      toType = toType,
      OrgDb = OrgDb
    )
    #row.names(eg)=eg$SYMBOL
    if (showConvertedGeneId) {
      cat("\n\n")
      print(knitr::kable(eg))
      cat("\n\n")
    }
    selectedProteinsToEg = unique(eg$ENTREZID)

    #If universe provided, universe may need ID transfer too.
    if (!is.null(universe)) {
      egUniverse = bitr(
        universe,
        fromType = fromType,
        toType = toType,
        OrgDb = OrgDb
      )
      universe=unique(egUniverse$ENTREZID)
    }
  } else {
    selectedProteinsToEg = selectedProteins
  }

  #browser()
  if ("Reactome" %in% modules) {
    if(showProgress){
      print("Running Reactome")
    }
    if (organism=="hsa") {
      organismReactome="human"
    } else if (organism=="mmu") {
      organismReactome="mouse"
    } else if (organism=="rno") {
      organismReactome="rat"
    } else {
      stop("organism Not supported for Reactome")
    }
    xReactome <- enrichPathway(
      gene = selectedProteinsToEg,
      organism     = organismReactome,
      universe = universe,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      readable = TRUE
    )
    xReactome = xReactome %>% 
      dplyr::arrange((pvalue)) %>% 
      dplyr::mutate(qscore = -log(qvalue, base = 10),
                    NegLog10pAdj = -log(p.adjust, base = 10),
                    CountRatio = sapply(strsplit(GeneRatio, "/"), function(x)
                      as.integer(x[1]) / as.integer(x[2]))
                   )
    listEnrichments[["Reactome"]] = xReactome
  }

  if ("KEGG" %in% modules) {
    if(showProgress){
      print("Running KEGG")
    }
    kk <- enrichKEGG(
      gene         = selectedProteinsToEg,
      organism     = organism,
      universe = universe,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff
    )
    kkx <- setReadable(kk, OrgDb, toType)
    kkx = kkx %>% 
      dplyr::arrange((pvalue)) %>% 
      dplyr::mutate(qscore = -log(qvalue, base = 10),
                    NegLog10pAdj = -log(p.adjust, base = 10),
                    CountRatio = sapply(strsplit(GeneRatio, "/"), function(x)
                      as.integer(x[1]) / as.integer(x[2]))
                    )
    listEnrichments[["KEGG"]] = kkx
  }

  if ("GO.BP" %in% modules) {
    if(showProgress){
      print("Running GO.BP")
    }
    ego <- enrichGO(
      gene          = selectedProteinsToEg,
      OrgDb         = OrgDb,
      universe = universe,
      ont           = "BP",
      pAdjustMethod = "BH",
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      readable      = TRUE
    )
    ego = ego %>% 
      dplyr::arrange((pvalue)) %>% 
      dplyr::mutate(qscore = -log(qvalue, base = 10),
                    NegLog10pAdj = -log(p.adjust, base = 10),
                    CountRatio = sapply(strsplit(GeneRatio, "/"), function(x)
                      as.integer(x[1]) / as.integer(x[2]))
                    )
    listEnrichments[["GO.BP"]] = ego
  }

  if ("GO.MF" %in% modules) {
    if(showProgress){
      print("Running GO.MF")
    }
    ego <- enrichGO(
      gene          = selectedProteinsToEg,
      OrgDb         = OrgDb,
      universe = universe,
      ont           = "MF",
      pAdjustMethod = "BH",
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      readable      = TRUE
    )
    ego = ego %>% 
      dplyr::arrange((pvalue)) %>% 
      dplyr::mutate(qscore = -log(qvalue, base = 10),
                    NegLog10pAdj = -log(p.adjust, base = 10),
                    CountRatio = sapply(strsplit(GeneRatio, "/"), function(x)
                      as.integer(x[1]) / as.integer(x[2]))
                    )
    listEnrichments[["GO.MF"]] = ego
  }
  if ("GO.CC" %in% modules) {
    if(showProgress){
      print("Running GO.CC")
    }
    ego <- enrichGO(
      gene          = selectedProteinsToEg,
      OrgDb         = OrgDb,
      universe = universe,
      ont           = "CC",
      pAdjustMethod = "BH",
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      readable      = TRUE
    )
    ego = ego %>% 
      dplyr::arrange((pvalue)) %>% 
      dplyr::mutate(qscore = -log(qvalue, base = 10),
                    NegLog10pAdj = -log(p.adjust, base = 10),
                    CountRatio = sapply(strsplit(GeneRatio, "/"), function(x)
                      as.integer(x[1]) / as.integer(x[2]))
                    )
    listEnrichments[["GO.CC"]] = ego
  }

  if ("WikiPathways" %in% modules) {
    if(showProgress){
      print("Running WikiPathways")
    }
    #get_wp_organisms()
    if (organism=="hsa") {
      organismWP="Homo sapiens"
    } else if (organism=="mmu") {
      organismWP="Mus musculus"
    } else if (organism=="rno") {
      organismWP="Rattus norvegicus"
    } else {
      stop("organism Not supported for WikiPathways")
    }
    wp <- enrichWP(
      gene          = selectedProteinsToEg,
      organism         = organismWP,
      universe = universe,
      pAdjustMethod = "BH",
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff
    )
    wpx <- setReadable(wp, OrgDb, toType) #transform geneID to input genesymbolID
    wpx = wpx %>% 
      dplyr::arrange((pvalue)) %>% 
      dplyr::mutate(qscore = -log(qvalue, base = 10),
                    NegLog10pAdj = -log(p.adjust, base = 10),
                    CountRatio = sapply(strsplit(GeneRatio, "/"), function(x)
                      as.integer(x[1]) / as.integer(x[2]))
                    )
    listEnrichments[["WikiPathways"]] = wpx
  }

  for(pathway in names(listEnrichments)){
    listEnrichments[[pathway]]@result=add_coverage(listEnrichments[[pathway]]@result)
  }
  return(listEnrichments)
}


makeDotPlotClusterProfilerEnrichment=function(dp,
                                              minCount=5,
                                              top=10,
                                              valueCut=NULL,
                                              y="Coverage",
                                              returnTable=FALSE,
                                              showProgress=TRUE,
                                              warp_len=50) {
  topBy="pvalue"
  colorBy="NegLog10pAdj"
  colorTitle="-log10(FDR)"

  if(showProgress){
    print(paste0("Pathways were filtered based on Count >= ",minCount, ifelse(is.null(valueCut), "", paste0(" and ", topBy, "<=", valueCut))))
  }
  if (class(dp)=="list") { #list as input, make into data.frame
    dataForPlotFrame=NULL
    for (i in 1:length(dp)) {
      temp=data.frame(dp[[i]],Module=names(dp)[i])
      dataForPlotFrame=rbind(dataForPlotFrame,temp)
    }
  } else { #data frame, extract needed columns
    if ("Module" %in% colnames(dp)) {
      dataForPlotFrame=dp[,unique(c("Description","Count","CountRatio",topBy,colorBy,y,"Module"))]
    } else {
      dataForPlotFrame=data.frame(dp[,unique(c("Description","Count","CountRatio",topBy,colorBy,y))],Module="1")
    }
  }
  #browser()
  temp=apply(dataForPlotFrame,2,function(x) all(is.na(x)))
  if (any(temp)) {
    warning(paste0("All NA values in ",names(which(temp)),". Skip making figure"))
    return(NULL)
  }

  if(!is.null(valueCut)){
    dataForPlotFrameP=dataForPlotFrame %>% 
      dplyr::group_by(Module) %>% 
      dplyr::filter(Count>=minCount & !!sym((topBy))<=valueCut)
  }else{
    dataForPlotFrameP=dataForPlotFrame %>% 
      dplyr::group_by(Module) %>% 
      dplyr::filter(Count>=minCount)
  }
  if (nrow(dataForPlotFrameP) == 0) {
    warning(paste0("No valid pathway. Skip making figure"))
    return(NULL)
  }
  
  dataForPlotFrameP = dataForPlotFrameP %>% 
    dplyr::group_by(Module) %>%
    dplyr::top_n(-top,wt=get(topBy)) %>%
    dplyr::arrange(desc(!!sym(y))) %>%
    as.data.frame()

  dataForPlotFrameP$DescriptionWrap=gsub("_+[a-z0-9 \\.A-Z\\(\\)-]+$","",dataForPlotFrameP$Description)
  dataForPlotFrameP$DescriptionWrap=str_wrap(dataForPlotFrameP$DescriptionWrap, width = warp_len)
  dataForPlotFrameP$DescriptionWrap=factor(dataForPlotFrameP$DescriptionWrap, levels=rev(unique(dataForPlotFrameP$DescriptionWrap)))

  max_value=max(dataForPlotFrameP[[colorBy]],na.rm=TRUE)

  p=dataForPlotFrameP %>%
    ggplot(aes(x=DescriptionWrap,y=Coverage,colour=!!sym(colorBy)))+
    geom_point(aes(size=Count)) +
    coord_flip() + 
    xlab("") +
    scale_colour_continuous(limits=c(0,max_value), low='blue', high='red',guide=guide_colorbar(reverse=FALSE))+
    guides(colour=guide_legend(title=colorTitle)) +
    facet_wrap(~Module,scales = "free")+
    theme_bw()+
    theme(
      #panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill="white")
    )
  if (returnTable) {
    return(list(
      figure=p,
      result=dataForPlotFrameP |> dplyr::select(-DescriptionWrap)
    ))
  } else {
    return(p)
  }
}

makeBarPlotClusterProfilerEnrichment=function(dp,minCount=5,top=10,valueCut=NULL,
                                              topBy="pvalue",fill="p.adjust",y="Coverage",
                                              returnTable=FALSE,
                                              showProgress=TRUE,
                                              warp_len=50) {
  if(showProgress){
    print(paste0("Pathways were filtered based on Count>=",minCount," and ",fill,"<=",ifelse(is.null(valueCut), "NULL", valueCut)))
  }
  if (class(dp)=="list") { #list as input, make into data.frame
    dataForPlotFrame=NULL
    for (i in 1:length(dp)) {
      temp=data.frame(dp[[i]],Module=names(dp)[i])
      dataForPlotFrame=rbind(dataForPlotFrame,temp)
    }
  } else { #data frame, extract needed columns
    if ("Module" %in% colnames(dp)) {
      dataForPlotFrame=dp[,unique(c("Description","Count","CountRatio",topBy,fill,y,"Module"))]
    } else {
      dataForPlotFrame=data.frame(dp[,unique(c("Description","Count","CountRatio",topBy,fill,y))],Module="1")
    }
  }
  #browser()
  temp=apply(dataForPlotFrame,2,function(x) all(is.na(x)))
  if (any(temp)) {
    warning(paste0("All NA values in ",names(which(temp)),". Skip making figure"))
    return(NULL)
  }

  if(!is.null(valueCut)){
    dataForPlotFrameP=dataForPlotFrame %>% 
      dplyr::group_by(Module) %>% 
      dplyr::filter(Count>=minCount & !!sym((fill))<=valueCut)
  }else{
    dataForPlotFrameP=dataForPlotFrame %>% 
      dplyr::group_by(Module) %>% 
      dplyr::filter(Count>=minCount)
  }
  if (nrow(dataForPlotFrameP) == 0) {
    warning(paste0("No valid pathway. Skip making figure"))
    return(NULL)
  }
  
  dataForPlotFrameP = dataForPlotFrameP %>% 
    dplyr::top_n(-top,wt=get(topBy)) %>%
    dplyr::mutate(Description = reorder_within(x=Description,by=!!sym(y), within=Module)) %>%
    dplyr::arrange(desc(!!sym(y)))

  dataForPlotFrameP$Description=gsub("_+[a-z0-9 \\.A-Z\\(\\)-]+$","",dataForPlotFrameP$Description)

  p=dataForPlotFrameP %>%
    ggplot(aes_string(x="Description",y=y,fill=fill))+geom_bar(stat = "identity")+
        scale_x_reordered()+
        scale_x_discrete(labels = function(x) str_wrap(x, width = warp_len))+
        coord_flip() + xlab("")
  if (min(dataForPlotFrame[,fill],na.rm=TRUE)<0) { #has negative values
    p=p+scale_fill_continuous(low='blue', high='red',guide=guide_colorbar(reverse=TRUE))
  } else {
    p=p+scale_fill_gradient(low='red', high='blue',guide=guide_colorbar(reverse=TRUE))
  }
  p=p+facet_wrap(~Module,scales = "free")+
    theme_bw()+
    theme(
      #panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  if (returnTable) {
    return(list(
      figure=p,
      result=dataForPlotFrameP
    ))
  } else {
    return(p)
  }
}


makeHeatmapClusterProfilerEnrichment=function(dataForPlot,
                                              minCount=5,
                                              top=3,
                                              topBy="pvalue",
                                              fill="NegLog10pAdj",
                                              returnData=FALSE,
                                              cluster_rows =TRUE,
                                              cluster_columns=TRUE,
                                              changeNA=NULL,
                                              na_col="grey") { #Please note changeNA values will change heatmap cluster

  if (class(dataForPlot)=="list") { #list as input, make into data.frame
    dataForPlotFrame=NULL
    for (i in 1:length(dataForPlot)) {
      if (!is.null(dataForPlot[[i]])) {
        temp=data.frame(dataForPlot[[i]],Module=names(dataForPlot)[i])
        dataForPlotFrame=rbind(dataForPlotFrame,temp)
      }
    }
  } else { #data frame, extract needed columns
    if (class(dataForPlot)=="compareClusterResult") {
      dataForPlotFrame=as.data.frame(dataForPlot[,c("Description","setSize",fill,topBy,"Cluster")])
      dataForPlotFrame=dataForPlotFrame %>% rename("setSize" = "Count",
                                                  "Cluster" = "Module")
    } else {
      dataForPlotFrame=dataForPlot[,c("Description","Count","CountRatio",topBy,"Module")]
    }
  }

  selectedPathways=dataForPlotFrame %>% 
    dplyr::group_by(Module) %>% 
    dplyr::filter(Count>=minCount) %>% 
    dplyr::top_n(-top,wt=get(topBy)) %>% 
    dplyr::pull(Description)

  dataForPlotFrameForHeatmap=dataForPlotFrame %>% 
    dplyr::filter(Description %in% selectedPathways) %>%
    dplyr::select("Description",all_of(fill),"Module") %>% 
    dplyr::pivot_wider(names_from = Module, values_from = all_of(fill))

  temp=colnames(dataForPlotFrameForHeatmap)
  dataForPlotFrameForHeatmap=data.frame(dataForPlotFrameForHeatmap)
  colnames(dataForPlotFrameForHeatmap)=temp #keep the colnames with different chracters

  row.names(dataForPlotFrameForHeatmap)=dataForPlotFrameForHeatmap[,1]
  dataForPlotFrameForHeatmap=as.matrix(dataForPlotFrameForHeatmap[,-1])
  if (!is.null(changeNA)) {
    dataForPlotFrameForHeatmap[is.na(dataForPlotFrameForHeatmap)]=changeNA
  }
  if (fill=="NegLog10pAdj") {
    heatmapTitle="-log10(pAdjusted)"
  } else {
    heatmapTitle=fill
  }
  if (returnData) {
    p=dataForPlotFrameForHeatmap
  } else {
    dataForPlotFrameForHeatmapMaxAbsValue=max(abs(dataForPlotFrameForHeatmap),na.rm=TRUE)
    dataForPlotFrameForHeatmapMinValue=min(dataForPlotFrameForHeatmap,na.rm=TRUE)
    if (dataForPlotFrameForHeatmapMinValue<0) { #has negative values
      col_fun = circlize::colorRamp2(c(-dataForPlotFrameForHeatmapMaxAbsValue, 0, dataForPlotFrameForHeatmapMaxAbsValue), c("green", "white", "red"))
    } else {
      col_fun = circlize::colorRamp2(c(0, dataForPlotFrameForHeatmapMinValue, dataForPlotFrameForHeatmapMaxAbsValue), c("white", "white", "red"))
    }
    p=Heatmap(dataForPlotFrameForHeatmap, 
              heatmap_legend_param = list(title = heatmapTitle),
              col=col_fun,
              na_col = na_col,
              cluster_rows =cluster_rows,
              cluster_columns=cluster_columns)
  }
  return(p)
}


#for a small number of selected proteins among for small number of universe proteins.
#For example, 30 proteins selected in 70 proteins universe.
makePercentInUniverseBarPlotClusterProfilerEnrichment=function(enrichmentResult,minCount=5,module="KEGG") {
  dataForPlot=data.frame(enrichmentResult[[module]]) %>%
    dplyr::filter(Count>=minCount) %>%
    dplyr::mutate(Selected=Count/sapply(strsplit(BgRatio,"/"),function(x) as.numeric(x[1])),
                  BgNotSelected=1-Selected) %>% 
    dplyr::mutate(Description=forcats::fct_reorder(Description, Selected))

  dataForPlot=dataForPlot[,c("Description","Selected","BgNotSelected")] %>%
    pivot_longer(cols =c("Selected","BgNotSelected"),names_to = "Category",values_to = "Proportion")

  p=ggplot(dataForPlot,aes(y=Description,x=Proportion ,fill=Category))+
    geom_col()+
    xlab("Proportion")+
    ggtitle(module)

  return(p)
}


enrichmentObjToTable=function(dataForPlot,splitGenes=FALSE) {
  if (class(dataForPlot)=="list") { #list as input, make into data.frame
    dataForPlotFrame=NULL
    for (i in 1:length(dataForPlot)) {
      if (!is.null(dataForPlot[[i]])) {
        temp=data.frame(dataForPlot[[i]],Module=names(dataForPlot)[i])
        dataForPlotFrame=rbind(dataForPlotFrame,temp)
      }
    }
  } else { #data frame, extract needed columns
    otherInfoColumns=intersect(colnames(head(dataForPlot)),c("pvalue","NegLog10pAdj","geneID"))

    if (class(dataForPlot)=="compareClusterResult") {
      dataForPlotFrame=as.data.frame(dataForPlot[,c("Description","setSize",otherInfoColumns,"Cluster")])
      dataForPlotFrame=dataForPlotFrame %>% rename("setSize" = "Count",
                                                   "Cluster" = "Module")
    } else {
      if ("Module" %in% colnames(dataForPlot)) {
        dataForPlotFrame=dataForPlot[,c("Description","Count","CountRatio",otherInfoColumns,"Module")]
      } else {
        #dataForPlotFrame=data.frame(dataForPlot[,c("Description","Count","CountRatio",otherInfoColumns)],Module="")
        dataForPlotFrame=data.frame(dataForPlot[,c("Description","Count","CountRatio",otherInfoColumns)])
      }
    }
  }

  if (splitGenes) {
    dataForPlotFrame <- separate_rows(dataForPlotFrame, "geneID", sep = "/")
  }
  return(dataForPlotFrame)
}


showGenesValueClusterProfilerEnrichment=function( dataForPlot,
                                                  geneValue,
                                                  minCount=5,
                                                  top=3,
                                                  selectedPathways=NULL,
                                                  topBy="pvalue",
                                                  fill="NegLog10pAdj",
                                                  returnData=FALSE,
                                                  rankGeneByModule=NULL,
                                                  proteinToOnePathway=FALSE,
                                                  changeNA=NULL,
                                                  na_col="grey",
                                                  figType=c("boxplot","dotplot")) { #Please note changeNA values will change heatmap cluster
  figType=match.arg(figType)
  dataForPlotFrame=enrichmentObjToTable(dataForPlot,splitGenes = TRUE)

  if (is.null(selectedPathways)) {
    #selectedPathways=unique(dataForPlotFrame %>% dplyr::select(-geneID)) %>% group_by(Module) %>% filter(Count>=minCount) %>% top_n(-top,wt=get(topBy)) %>% pull(Description)
    selectedPathways=unique(dataForPlotFrame %>% 
      dplyr::select(-geneID)) %>% 
      dplyr::filter(Count>=minCount) %>% 
      dplyr::top_n(-top,wt=get(topBy)) %>% 
      dplyr::pull(Description)

    selectedPathways=unique(selectedPathways)
  }



  if (proteinToOnePathway) {
    dataForPlotFrame=dataForPlotFrame %>% 
      dplyr::group_by(geneID) %>% 
      dplyr::slice_min(pvalue) %>% 
      dplyr::ungroup()
  }

  if ("Module" %in% colnames(geneValue) & "Module" %in% colnames(dataForPlotFrame)) {
    dataForPlotFrameWithValue=dataForPlotFrame %>% 
      dplyr::left_join(geneValue,by=c("geneID","Module"))
  } else {
    #dataForPlotFrameWithValue=dataForPlotFrame %>% left_join(geneValue,by=c("geneID"))
    dataForPlotFrameWithValue=dataForPlotFrame %>% 
      dplyr::inner_join(geneValue,by=c("geneID"))
  }


  dataForPlotFrameWithValue=dataForPlotFrameWithValue %>% 
    dplyr::filter(Description %in% selectedPathways)

  valueColumnName=setdiff(colnames(geneValue),c("geneID","Module"))[1]

  #browser()
  DescriptionValues=tapply(dataForPlotFrameWithValue[[valueColumnName]],dataForPlotFrameWithValue$Description,median,na.rm=TRUE)
  dataForPlotFrameWithValue$Description <- factor(dataForPlotFrameWithValue$Description, levels =names(sort(DescriptionValues)) )

  #browser()
  #valueColumnName=enquo(valueColumnName)
  if (figType=="boxplot") {
    p=ggplot(dataForPlotFrameWithValue,aes(x=!!sym(valueColumnName),y=Description))+geom_boxplot()
  } else if (figType=="dotplot") {

    #make geneID column as a factor and levels ordered by median value of valueColumnName
    if (!is.null(rankGeneByModule)) {
      temp=dataForPlotFrameWithValue[which(dataForPlotFrameWithValue$Module==rankGeneByModule),]
    } else {
      temp=dataForPlotFrameWithValue
    }
    dataForPlotFrameWithValue$geneID <- factor(dataForPlotFrameWithValue$geneID, levels =names(sort(tapply(temp[[valueColumnName]],temp$geneID,median,na.rm=TRUE))) )
    p=ggplot(dataForPlotFrameWithValue,aes(x=geneID,y=!!sym(valueColumnName),color=Module))+geom_point()+facet_grid(.~Description,scale="free",space="free")+ylim(range(dataForPlotFrameWithValue[[valueColumnName]]))
  }
  return(p+theme_bw())
}

simplifyClusterProfilerEnrichment=function(enrichmentObj,method=c("pvalue","Count","CountRatio")) {
  method=match.arg(method)
  dataForPlotFrame=enrichmentObjToTable(enrichmentObj,splitGenes = TRUE)

  if (method=="pvalue") {
    dataForPlotFrameSimplified=dataForPlotFrame %>% 
      dplyr::group_by(geneID) %>% 
      dplyr::slice_min(pvalue)
  } else if (method=="Count" | method=="CountRatio") {
    dataForPlotFrameSimplified=dataForPlotFrame %>% 
      dplyr::group_by(geneID) %>% 
      dplyr::slice_max(Count)
  }
  #return(dataForPlotFrameSimplified[,c("geneID","Description")])
  return(structure(dataForPlotFrameSimplified[["Description"]],names=dataForPlotFrameSimplified[["geneID"]]))
}

compareClusterByClusterProfiler=function( geneList,
                                          fun=GSEA,
                                          idType=c("gene_symbol","entrez_gene","ensembl_gene"),
                                          msigdbrCat="C2",
                                          msigdbrSubCat=NULL,
                                          merge=FALSE) {
  require(msigdbr)
  require(DOSE)
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(enrichplot)

  idType=match.arg(idType)
  #print(idType)
  if (is.null(msigdbrSubCat)) {
    geneSetsToIdTable <- msigdbr(species = "Homo sapiens", category = msigdbrCat)
  } else {
    geneSetsToIdTable <- msigdbr(species = "Homo sapiens", category = msigdbrCat,subcategory=msigdbrSubCat)
  }
  geneSetsToIdTable=geneSetsToIdTable%>% dplyr::select(gs_name, all_of(idType))
  #print(head(geneSetsToIdTable))

  ##############################
  # use loop to do enrichment, so that individual results will be saved
  ##############################
  #geneListInSetsResult <- compareCluster(geneCluster = geneList, fun = fun, TERM2GENE = geneSetsToIdTable)
  # if (idType=="entrez_gene") {
  #   geneListInSetsResult <- setReadable(geneListInSetsResult, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  # }
  geneListInSetsResult=list()
  for (i in 1:length(geneList)) {
    geneListInSetsResult[[i]]=fun(geneList[[i]],TERM2GENE = geneSetsToIdTable,pvalueCutoff = 1)
  }
  names(geneListInSetsResult)=names(geneList)

  if (merge) {
    geneListInSetsResult=merge_result(geneListInSetsResult)
  }

  return(geneListInSetsResult)
}

compareClusterVis=function(geneListInSetsResult,plotList=c("dotplot","cnetplot","emapplot")) {
  p1=dotplot(geneListInSetsResult)
  p2=cnetplot(geneListInSetsResult,pie.params = list(pie = "count"))

  temp=pairwise_termsim(geneListInSetsResult)
  p3=emapplot(temp,pie.params = list(pie = "count"))
  p4=treeplot(temp)

  #heatplot(geneListInSetsResult, showCategory=5)
  #upsetplot(geneListInSetsResult, showCategory=5)

  #barplot(geneListInSetsResult)

}

do_pathway<-function( genes, 
                      universe = NULL,
                      modules = c("Reactome","KEGG"),
                      organism = "hsa",
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = "org.Hs.eg.db") {
  dataForPlot=enrichmentByClusterProfiler(
    unique(genes),
    universe=universe,  #universe=unique(ls_res$all$gene),
    modules=modules,
    organism=organism,
    fromType=fromType,
    toType=toType,
    OrgDb=OrgDb,
    minGSSize = 10,
    maxGSSize = 1000,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    showConvertedGeneId=FALSE,
    showProgress=FALSE)

  for(gname in names(dataForPlot)){
    dataForPlot[[gname]]@result$Description=gsub('\\(NMD\\)|\\(EJC\\)','',dataForPlot[[gname]]@result$Description)
  }

  num_background_genes = ifelse(all(is.null(universe)), 'all', length(unique(universe)))
  return(list(dataForPlot=dataForPlot, num_genes=length(unique(genes)), num_background_genes=num_background_genes))
}

get_pathway_table<-function(pathway_res, csv_file=NULL,qvalueCut=0.05){
  res_tbl = pathway_res[pathway_res$qvalue <qvalueCut ,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "qvalue", "Count", "CountRatio", "Coverage")]
  if(!is.null(csv_file)){
    write.csv(res_tbl, csv_file, row.names=F)
  }
  return(list(res_tbl=res_tbl, csv_file=csv_file))
}

get_pathway_figure<-function(dp, png_file, fig_width, fig_height, pname=NULL, y="Coverage", ylab="Geneset coverage (%)", plot_type="bar"){
  if(!is.null(pname)){
    dp = list(dp[[pname]])
    names(dp) = pname
  }
  if(plot_type == "dot"){
    g = makeDotPlotClusterProfilerEnrichment(dp, top=20, valueCut=NULL, y=y, showProgress=FALSE)
  }else{
    g = makeBarPlotClusterProfilerEnrichment(dp, top=20, valueCut=NULL, y=y, showProgress=FALSE)
  }
  if(is.null(g)){
    return(NULL)
  }
  
  g=g + ylab(ylab) +
    theme(text = element_text(size=8), strip.background =element_rect(fill="white"), strip.text = element_text(size=8, face="bold"))

  ggsave(png_file, g, width=fig_width, height=fig_height, units="in", dpi=300, bg="white")
  return(list(g=g, png_file=png_file))
}

show_pathway<-function( dataForPlotList, 
                        pname, 
                        prefix="pathway",
                        qvalueCut=0.05,
                        pathway_width=6, 
                        pathway_height=6, 
                        y="Coverage", 
                        ylab="Geneset coverage (%)",
                        plot_type="bar"){
  prefix_name = paste0(prefix, ".", pname)

  tbl=get_pathway_table(
    pathway_res=dataForPlotList[[pname]]@result,
    qvalueCut=qvalueCut,
    csv_file=paste0(prefix_name, ".csv"))

  if(nrow(tbl$res_tbl) > 0) {
    dp_new = dataForPlotList[[pname]]
    dp_new@result = dp_new@result |> dplyr::filter(qvalue < qvalueCut)
    dp = list(dp_new)
    names(dp) = pname    
  }else{
    dp=dataForPlotList
  }

  fig_height=ifelse(nrow(tbl$res_tbl) >= 15, pathway_height * 1.5, pathway_height)
  
  fig=get_pathway_figure(
    dp=dp, 
    png_file=paste0(prefix_name, ".png"), 
    fig_width=pathway_width, 
    fig_height=fig_height, 
    pname=pname,
    y=y, 
    ylab=ylab,
    plot_type=plot_type)
  
  if(is.null(fig)){
    return(list(res_tbl=tbl$res_tbl, png=NULL, csv=tbl$csv_file))
  }else{
    return(list(res_tbl=tbl$res_tbl, png=fig$png_file, csv=tbl$csv_file))
  }
}

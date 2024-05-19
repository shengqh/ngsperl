library(tibble)

add_secreted<-function(df){
  selectedProteinsCategory=rep("Intracellular", nrow(df))
  selectedProteinsCategory[grep("membrane proteins", df[["Protein.class"]])] = "Membrane"
  selectedProteinsCategory[grep("secreted proteins", df[["Protein.class"]])] = "Secreted"
  df$Secreted.Category = selectedProteinsCategory
  return(df)
}

findProteinInDatabase = function(proteinDatabaseAll, selectedProteins,selectedProteinsLabels=NULL) {
  proteinDatabaseAllOtherNamesList = strsplit(proteinDatabaseAll$`Gene synonym`, ",\ ")
  
  selectedProteinsMatch = data.frame()
  for (i in seq_along(selectedProteins)) {
    selectedProteinOne = selectedProteins[i]
    selectedProteinOneInDbInd = which(proteinDatabaseAll$Gene == selectedProteinOne)
    if (length(selectedProteinOneInDbInd) > 0) {
      #matching "Gene" column
      selectedProteinsMatch = rbind(
        selectedProteinsMatch,
        c(
          selectedProteins[i],
          selectedProteinOneInDbInd[1],
          proteinDatabaseAll[selectedProteinOneInDbInd[1], ][[1]]
        )
      )
    } else {
      #NOT matching "Gene" column
      selectedProteinOneInDbInd = which(
        sapply(proteinDatabaseAllOtherNamesList, function(x)
          selectedProteinOne %in% x)
      )
      if (length(selectedProteinOneInDbInd) > 0) {
        selectedProteinsMatch = rbind(
          selectedProteinsMatch,
          c(
            selectedProteins[i],
            selectedProteinOneInDbInd[1],
            proteinDatabaseAll[selectedProteinOneInDbInd[1], ][[1]]
          )
        )
      } else {
        selectedProteinsMatch = rbind(selectedProteinsMatch, c(selectedProteins[i], NA, NA))
      }
    }
  }
  colnames(selectedProteinsMatch) = c("GeneSymbol", "MatchedInd", "GeneNameInDb")
  if (!is.null(selectedProteinsLabels)) {
    selectedProteinsMatch$GeneLabelForFigure = selectedProteinsLabels
  } else {
    selectedProteinsMatch$GeneLabelForFigure = selectedProteinsMatch$GeneSymbol
  }
  
  proteinDatabase = proteinDatabaseAll[selectedProteinsMatch[, "MatchedInd"], ]
  
  print(paste0(length(selectedProteins)," proteins used in mapping"))
  selectedProteinsNotInGeneColumn = selectedProteinsMatch[which(is.na(selectedProteinsMatch[, 2])), 1]
  if (length(selectedProteinsNotInGeneColumn) > 0) {
    warning(paste0(length(selectedProteinsNotInGeneColumn)," proteins can't be mapped to Protein expression database"))
    print("Proteins NOT in Protein expression database:")
    print(selectedProteinsNotInGeneColumn)
  }
  
  return(list(proteinDatabase,selectedProteinsMatch))
  
}

findProteinInDatabase_fast = function(proteinDatabaseAll, selectedProteins, selectedProteinsLabels=NULL, showUnmappedGenes=FALSE) {
  proteinDatabaseAllOtherNamesList = strsplit(proteinDatabaseAll$`Gene synonym`, ",\ ")

  found_genes = intersect(selectedProteins, proteinDatabaseAll$Gene)
  found_genes_df = proteinDatabaseAll[which(proteinDatabaseAll$Gene %in% found_genes),]
  found_genes_df = found_genes_df[!duplicated(found_genes_df$Gene),]
  found_genes_df = found_genes_df %>% add_column(GeneSymbol = found_genes_df$Gene, .before="Gene")
  
  notfound_genes = setdiff(selectedProteins, proteinDatabaseAll$Gene)
  if(length(notfound_genes) > 0) {
    synonyms_genes = intersect(notfound_genes, unlist(proteinDatabaseAllOtherNamesList))
    if(length(synonyms_genes) > 0) {
      synonyms_genes_idx = unlist(lapply(synonyms_genes, function(x) {
        which(sapply(proteinDatabaseAllOtherNamesList, function(y) x %in% y))[1]}))
      synonyms_genes_df = proteinDatabaseAll[synonyms_genes_idx,,drop=FALSE]
      synonyms_genes_df = synonyms_genes_df %>% add_column(GeneSymbol = synonyms_genes, .before="Gene")
      found_genes_df = rbind(found_genes_df, synonyms_genes_df)
    }
  }

  final_df=merge(data.frame("GeneSymbol"=selectedProteins), found_genes_df, by="GeneSymbol", all.x=TRUE, fill=NA)
  if (!is.null(selectedProteinsLabels)) {
    final_df$GeneLabelForFigure = selectedProteinsLabels
  } else {
    final_df$GeneLabelForFigure = final_df$GeneSymbol
  }
  final_df=final_df[!is.na(final_df$Gene),]

  selectedProteinsNotInGeneColumn = setdiff(selectedProteins, found_genes_df$GeneSymbol)
  if (length(selectedProteinsNotInGeneColumn) > 0) {
    warning(paste0(length(selectedProteinsNotInGeneColumn)," out of ", length(selectedProteins), " genes can't be mapped to Protein expression database"))
    if(showUnmappedGenes){
      print("Proteins NOT in Protein expression database:")
      print(selectedProteinsNotInGeneColumn)      
    }
  }

  colnames(final_df) = make.names(colnames(final_df))
  final_df = add_secreted(final_df)

  return(list(final_df, selectedProteinsNotInGeneColumn))
}

summarizeCategoricalData<-function(df, category, maxCategory = 10, columnName = "RNA.tissue.specificity"){
  #use variable name as column name in group_by
  dataSelected = df %>% 
    dplyr::group_by(!!sym(columnName)) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::mutate(Category = category) %>% 
    dplyr::arrange(desc(n))
  dataSelected[which(is.na(dataSelected[, 1])), 1] = "Not Available"
  
  #keep top 10 rows (by n) in dataSelected
  if (nrow(dataSelected) > (maxCategory + 1)) {
    dataSelectedSub = rbind(dataSelected[1:maxCategory, ], c(
      "Others",
      n = sum(dataSelected[(maxCategory + 1):nrow(dataSelected), "n"]),
      Category = category
    ))
  } else {
    dataSelectedSub = dataSelected
  }
  dataSelectedSub$n = as.numeric(dataSelectedSub$n)
  return(list(dataSelected, dataSelectedSub))
}

summarizeCategoricalDataFromProteinDatabase = function( proteinDatabaseSelectedProteins,
                                                        proteinDatabaseUniverse,
                                                        maxCategory = 10,
                                                        columnName = "RNA.tissue.specificity") {
  #use variable name as column name in group_by
  temp = summarizeCategoricalData(
    df=proteinDatabaseSelectedProteins, 
    category="Selected", 
    maxCategory, 
    columnName)
  dataSelected = temp[[1]]
  dataSelectedSub = temp[[2]]
  rm(temp)
  
  if (!is.null(proteinDatabaseUniverse)) {
    temp = summarizeCategoricalData(proteinDatabaseUniverse, "Background", 1000, columnName)
    dataUniverse = temp[[1]]
    rm(temp)

    if (nrow(dataUniverse) > (maxCategory + 1)) {
      #we need to keep the category identical to the selected data
      d1 = dataUniverse[as.character(unlist(dataUniverse[, 1])) %in% as.character(unlist(dataSelectedSub[, 1])),,drop=FALSE ]
      d2 = dataUniverse[!(as.character(unlist(dataUniverse[, 1])) %in% as.character(unlist(dataSelectedSub[, 1]))),,drop=FALSE ]
      dataUniverseSub = rbind(d1, c(
        "Others",
        n = sum(d2[, "n"]),
        Category = "Background"
      ))
    } else {
      dataUniverseSub = dataUniverse
    }
    dataUniverseSub$n = as.numeric(dataUniverseSub$n)

    dataForPlot = rbind(dataSelectedSub, dataUniverseSub)
    dataForPlot$n = as.numeric(dataForPlot$n)
    
    dataForAnalysis = dataSelected %>% 
      left_join(dataUniverse, by = columnName) %>% 
      mutate(Enrichment = n.x / n.y) %>% 
      arrange(desc(Enrichment))

    nSelected = nrow(proteinDatabaseSelectedProteins)
    nUniverse = nrow(proteinDatabaseUniverse)
    
    #run a Hypergeometric test by phyper
    dataForAnalysis = dataForAnalysis %>% 
      mutate( pvalue = phyper(n.x - 1, nSelected, nUniverse - nSelected, n.y, lower.tail = FALSE),
              pAdj = p.adjust(pvalue, method = "BH"),
              pAdjScore = -log10(pAdj))
    
    dataForAnalysis=dataForAnalysis %>% 
      slice_max(n=maxCategory,order_by=pAdjScore) %>% 
      mutate(Levels=factor(!!sym(columnName),levels=unique(!!sym(columnName))[order(pAdjScore)]))

    return(list(dataForPlot, dataForAnalysis))
  } else {
    return(list(dataSelectedSub))
  }
}

summarizenTPMDataFromProteinDatabase = function(proteinDatabaseSelectedProteins,
                                                proteinDatabaseUniverse,
                                                maxCategory = 20,
                                                columnName = "RNA.tissue.specific.nTPM") {
  columnTypeName = gsub("specific.nTPM", "specificity", columnName)
  
  #split data into multiple rows
  TissueEnrichSelectedProteins = proteinDatabaseSelectedProteins[, c("Gene", columnTypeName, columnName)] %>% separate_rows(!!sym(columnName),
                                                                                                                            sep = ";") %>% separate(!!sym(columnName),
                                                                                                                                                    into = c("Tissue", "nTPM"),
                                                                                                                                                    sep = ": ")
  TissueEnrichSelectedProteins = na.omit(TissueEnrichSelectedProteins)
  
  print(knitr::kable(unique(TissueEnrichSelectedProteins[,c("Gene",columnTypeName)]) %>% group_by(!!sym(columnTypeName)) %>% summarise(Count=n()),caption = "Number of Proteins with Enriched/Enhanced expression"))
  
  tissueOrder = names(sort(table(TissueEnrichSelectedProteins$Tissue), decreasing = TRUE))
  
  TissueEnrichSelectedProteinsTypeMatrix = pivot_wider(
    TissueEnrichSelectedProteins[, -4],
    names_from = Tissue,
    values_from = !!sym(columnTypeName)
  )
  TissueEnrichSelectedProteinsTypeMatrix = TissueEnrichSelectedProteinsTypeMatrix[, tissueOrder]
  #order row of dataForPlot by values from first column to last column
  TissueEnrichSelectedProteinsTypeMatrix = TissueEnrichSelectedProteinsTypeMatrix[do.call(order, TissueEnrichSelectedProteinsTypeMatrix),]
  
  TissueEnrichSelectedProteinsnTPMMatrix = pivot_wider(
    TissueEnrichSelectedProteins[, -2],
    names_from = Tissue,
    values_from = nTPM,
    values_fn = as.numeric
  )
  TissueEnrichSelectedProteinsnTPMMatrix = TissueEnrichSelectedProteinsnTPMMatrix[, tissueOrder]
  #order row of dataForPlot by values from first column to last column
  TissueEnrichSelectedProteinsnTPMMatrix = TissueEnrichSelectedProteinsnTPMMatrix[do.call(order,-TissueEnrichSelectedProteinsnTPMMatrix),]
  
  #ComplexHeatmap::Heatmap(scale(log(TissueEnrichSelectedProteinsnTPMMatrix)),cluster_rows = FALSE,cluster_columns = FALSE)
  
  if (exists("proteinDatabaseUniverse")) {
    #split data into multiple rows
    TissueEnrichUniverse = proteinDatabaseUniverse[, c("Gene", columnTypeName, columnName)] %>% separate_rows(!!sym(columnName),
                                                                                                              sep = ";") %>% separate(!!sym(columnName),
                                                                                                                                      into = c("Tissue", "nTPM"),
                                                                                                                                      sep = ": ")
    TissueEnrichUniverse = na.omit(TissueEnrichUniverse)
    TissueEnrichUniverseTypeMatrix = pivot_wider(
      TissueEnrichUniverse[, -4],
      names_from = Tissue,
      values_from = !!sym(columnTypeName)
    )
    
    temp1 = apply(TissueEnrichSelectedProteinsTypeMatrix, 2, function(x)
      length(which(!is.na(x))))
    temp2 = apply(TissueEnrichUniverseTypeMatrix[, colnames(TissueEnrichSelectedProteinsTypeMatrix)], 2, function(x)
      length(which(!is.na(x))))
    
    dataForAnalysis = data.frame(n.x = temp1, n.y = temp2) %>% rownames_to_column(var = columnName) %>%
      mutate(Enrichment = n.x / n.y) %>% arrange(desc(Enrichment)) %>% mutate(
        pvalue = phyper(n.x - 1, nSelected, nUniverse - nSelected, n.y, lower.tail =
                          FALSE),
        pAdj = p.adjust(pvalue, method = "BH"),
        pAdjScore = -log10(pAdj)
      )
    dataForAnalysis = dataForAnalysis %>% slice_max(n = maxCategory, order_by =
                                                      pAdjScore) %>% mutate(Levels = factor(!!sym(columnName), levels = unique(!!sym(columnName))[order(pAdjScore)]))
    return(
      list(
        TissueEnrichSelectedProteinsTypeMatrix,
        TissueEnrichSelectedProteinsnTPMMatrix,
        dataForAnalysis
      )
    )
  } else {
    return(
      list(
        TissueEnrichSelectedProteinsTypeMatrix,
        TissueEnrichSelectedProteinsnTPMMatrix
      )
    )
  }
}



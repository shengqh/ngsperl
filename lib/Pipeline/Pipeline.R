library(stringr)
library(ggplot2)
#library(tools)

display_webgestalt=function(files) {
  cat(paste0("  \n## WebGestalt  \n"))
  enrich_files<-files[grepl("WebGestalt_", rownames(files)),]
  enrich_files$Comparisons<-gsub("WebGestalt_GO....", "", enrich_files$V2)
  enrich_files$Comparisons<-gsub("WebGestalt_KEGG_", "", enrich_files$Comparisons)
  comparisons<-unique(enrich_files$Comparisons)
  for (j in 1:length(comparisons)){
    comparison<-comparisons[j]
    comp_files<-enrich_files[enrich_files$Comparisons == comparison,]
    cat(paste0("  \n### ", comparison, "  \n"))
  
    for (i in 1:nrow(comp_files)){
      if (!file.exists(comp_files[i,1])) {
        next;
      }
      ename<-rownames(comp_files)[i]
      if(grepl(".html.rds", comp_files[i,1])){
        plotData <-readRDS(comp_files[i,1])
        fdata<-plotData$enriched
        fdata<-fdata[c(1:min(nrow(fdata),5)),c(1,2,4,5,6,7,8,9,12,13)]
        print(kable(fdata, caption=tabRef(ename, ename)) %>% 
          kable_styling() %>%
          row_spec(which(fdata$geneUp > fdata$geneDown), color = "black", background = "bisque") %>% 
          row_spec(which(fdata$geneUp < fdata$geneDown), color = "black", background = "honeydew") %>%
          htmltools::HTML())
        cat("\n\n<hr>")
      }else{
        #txtFile<-gsub(".html.rds$", "", comp_files[i,1])
        #fdata<-read.table(txtFile, sep="\t", header=T, stringsAsFactors = F)
        fdata<-read.table(comp_files[i,1], sep="\t", header=T, stringsAsFactors = F)
        fdata<-fdata[c(1:min(nrow(fdata))),]
        fdata$geneSet<-paste0("[", fdata$geneSet, "](", fdata$link, ")")
        fdata<-fdata[,c(1,2,4,5,6,7,8,9)]
        print(kable(fdata, caption=tabRef(ename, ename)) %>%
        kable_styling() %>%
        htmltools::HTML())
      }
    }  
  }
}

processGseaTable=function(gseaTableFile,maxCategoryFdr=0.05,maxCategory=5,absLinkPath=FALSE) {
  rawTable<-read.delim(gseaTableFile,header=T,as.is=T)
  rawTableOut<-NULL
  for (i in head(which(rawTable$FDR.q.val<=maxCategoryFdr),maxCategory)) {
    categoryName=rawTable[i,1]
    gseaUrl=paste0("http://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=",categoryName)
    categoryNameInTable<-addLinkTag(text=categoryName,link=gseaUrl)
    
    gseaWeb<-getURL(gseaUrl)
    temp<-strsplit(gseaWeb,"description|\\<td\\>")[[1]]
    j=grep("Brief",temp)
    categoryDescription<-gsub("^>","",temp[j+2])
    categoryDescription<-gsub("<\\/$","",categoryDescription)
    
    reg_pattern = "<a href=.+?]</a>"

    if(length(categoryDescription) > 0){
      while(grepl(reg_pattern, categoryDescription)){
        url = str_extract(categoryDescription, reg_pattern) 
        
        name=str_match(url, '\\[(.+?)]')[[2]]
        name_url=str_match(url, "<a href='(.+)'")[[2]]
        namelink=addLinkTag(name, name_url)
        
        categoryDescription<-gsub(url, namelink, fixed=TRUE, categoryDescription)
      }
    }

    rawTableOut<-rbind(rawTableOut,unlist(c(categoryNameInTable,categoryDescription,rawTable[i,c(4,5,6,8)])))
  }
  if(!is.null(rawTableOut)){
    colnames(rawTableOut)[1:2]<-c("Name","Description")
  }
  return(rawTableOut)
}

print_gsea_subfolder=function(resultDirSub, gname, gseaCategory, is_positive, maxCategory=5, print_rmd=TRUE){
  pstr = ifelse(is_positive, "pos", "neg")
  title = ifelse(is_positive, "Positive-regulated", "Negative-regulated")
  pattern = paste0("gsea_report_for_na_", pstr, '_\\d+\\.(tsv|xls)$')
        
  gseaTableFile<-list.files(resultDirSub,pattern=pattern,full.names=TRUE)

  if (length(gseaTableFile)!=1) {
    warning(paste("Can't find ", title, "GSEA table file in", resultDirSub))
    return(NULL)
  }
  
  rawTableOut<-processGseaTable(gseaTableFile, maxCategory=maxCategory)
  if(is.null(rawTableOut)){
    warning(paste("Can't find significant", title, "GSEA gene set in", resultDirSub))
    return(NULL)
  }

  rawTable<-read.delim(gseaTableFile,header=T,as.is=T)
  rawTable<-rawTable[rawTable$FDR.q.val <= 0.05, c("NAME", "NES", "FDR.q.val")]

  if (nrow(rawTable) == 0){
    return(NULL)
  }

  rawTable$Core.Enrichment=0
  for(idx in c(1:nrow(rawTable))){
    gsfile = paste0(resultDirSub, "/", rawTable$NAME[idx], ".tsv")
    if(!file.exists(gsfile)){
      break
    }
    gsdetail = read.table(gsfile, sep="\t", header=T)
    rawTable$Core.Enrichment[idx]=sum(gsdetail$CORE.ENRICHMENT == "Yes")
  }
  return(list(rawTable=rawTable, rawTableOut=rawTableOut))
}

display_gsea=function(files, target_folder="", gsea_prefix="#", print_rmd=TRUE) {
  if(print_rmd){
    cat(paste0("\n\n", gsea_prefix, "# GSEA\n\n"))
  }

  result = data.frame(comparison=character(),
                      data_name=character(),
                      category=character(),
                      file_key=character(),
                      file_path=character())
  
  is_singlecell<-"compName" %in% colnames(files)

  maxCategory=ifelse(is_singlecell, 10, 5)
  
  if(is_singlecell){
    gsea_files<-files
    gsea_files$Comparisons<-gsea_files$compName
    comparisons<-rev(unique(gsea_files$Comparisons))
  }else{
    gsea_files<-files
    gsea_files$Comparisons<-gsea_files$V2
    comparisons<-unique(gsea_files$Comparisons)
  }
  j=2
  for (j in 1:length(comparisons)){
    comparison<-comparisons[j]
    comp_files<-gsea_files[gsea_files$Comparisons == comparison,]

    if(print_rmd){
      cat(paste0("\n\n", gsea_prefix, "## ", comparison, "\n\n"))
    }
  
    i=1
    for (i in 1:nrow(comp_files)){
      gname=comparison
      if(is_singlecell){
        gfolders<-comp_files
      }else{
        gfolders<-read.csv(comp_files[i,1],header=T,stringsAsFactor=F, check.names=F)
      }

      k=1
      for (k in 1:nrow(gfolders)){
        resultDirSub<-gfolders$Folder[k]
        gseaCategory<-gfolders$GseaCategory[k]
        pos = print_gsea_subfolder(resultDirSub, gname, gseaCategory, TRUE, maxCategory=maxCategory)
        neg = print_gsea_subfolder(resultDirSub, gname, gseaCategory, FALSE, maxCategory=maxCategory)
        if (is.null(pos) & is.null(neg)){
          next
        }

        prefix = paste0(gname, ".", gseaCategory)

        final = NULL
        if(!is.null(pos)){
          pos_file=paste0(target_folder, prefix, ".pos.csv")
          write.csv(pos$rawTableOut, pos_file)
          result[nrow(result) + 1,] = c(comparison, gname, gseaCategory, "pos_file", pos_file)

          final = pos$rawTable[pos$rawTable$Core.Enrichment > 0,]
        }
        if(!is.null(neg)){
          neg_file=paste0(target_folder, prefix, ".neg.csv")
          write.csv(neg$rawTableOut, neg_file)
          result[nrow(result) + 1,] = c(comparison, gname, gseaCategory, "neg_file", neg_file)

          final = rbind(final, neg$rawTable[neg$rawTable$Core.Enrichment > 0,])
        }

        if(is.null(final)){
          next
        }

        has_enriched = TRUE

        enriched_file=paste0(target_folder, prefix, ".enriched.csv")
        write.csv(final, enriched_file, row.names=F)
        result[nrow(result) + 1,] = c(comparison, gname, gseaCategory, "enriched_file", enriched_file)

        enriched_png=paste0(enriched_file, ".png")
        result[nrow(result) + 1,] = c(comparison, gname, gseaCategory, "enriched_png", enriched_png)

        final<-final[order(final$NES, decreasing=T),]
        final$NAME<-factor(final$NAME, levels=final$NAME)
        g<-ggplot(final, aes(NES, NAME)) + geom_point(aes(size=Core.Enrichment, color=FDR.q.val)) + geom_vline(xintercept=0) + theme_bw() + ylab("") + xlab("Normalized enrichment score")

        height=max(1000, 3000/40*nrow(final))
        warning(final$NAME)
        ncharmax=max(nchar(as.character(final$NAME)))
        width=max(3000, ncharmax * 40 + 1000)
        png(enriched_png, width=width, height=height, res=300)
        print(g)
        dev.off()

        if(print_rmd){
          cat(paste0("\n\n<img src='", pngfile, "'>\n\n"))
        }
        #include_graphics(file_path_as_absolute(pngfile))
        #cat("\n\n")
      }
    }  

    if(!has_enriched){
      if(print_rmd){
      	cat(paste0('\n#### <span style="color:red">WARNING: no geneset is significantly enriched.</span>\n\n'))
      }
    }
  }
  return(result)
}

get_versions=function(){
  files<-read.delim("fileList4.txt", header=F, stringsAsFactors = F)
  df<-NULL
  curfile<-files$V1[1]
  for(curfile in files$V1){
    if(file.exists(curfile)){
      curdf<-read.csv(curfile, header=F, check.names=F)
      df<-rbind(df, curdf)
    }
  }

  if(file.exists("fileList5.txt")){
    vers<-read.delim("fileList5.txt", header=F, stringsAsFactors = F)
    vers<-vers[,c(2,1)]
    colnames(vers)<-c("V1","V2")
    df<-rbind(df, vers)
  }

  df<-unique(df)
  df<-df[order(df$V1),]
  colnames(df)<-c("Software", "Version")
  return(df)
}

display_versions=function(){
  df<-get_versions()
  print(kable(df, caption=tabRef("versionFiles", "Software versions"), row.names=F) %>%
          kable_styling() %>%
          htmltools::HTML())
  return(df)
}

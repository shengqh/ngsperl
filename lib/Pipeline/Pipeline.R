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

    have_result = FALSE
    for (i in 1:nrow(comp_files)){
      if (!file.exists(comp_files[i,1])) {
        next;
      }
      have_result = TRUE
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
        fdata$geneSet<-paste0("[", fdata$geneSet, "](", fdata$link, "){target='_blank'}")
        fdata<-fdata[,c(1,2,4,5,6,7,8,9)]
        print(kable(fdata, caption=tabRef(ename, ename)) %>%
        kable_styling() %>%
        htmltools::HTML())
      }
    }

    if(!have_result){
      cat("\n\nNo WebGestalt result. It might caused by very limited number of differential expressed genes.\n<hr>")
    }
  }
}

processGseaTable=function(gseaTableFile, maxCategoryFdr=0.05) {
  if(require('msigdbr')){
    gs<-msigdbr:::msigdbr_genesets
    gsd_map=split(gs$gs_description, gs$gs_name)
  }else{
    gsd_map=list()
  }

  rawTable<-read.delim(gseaTableFile,header=T,as.is=T)
  rawTable<-subset(rawTable, NES != "---")

  rawTable<-rawTable[rawTable$FDR.q.val<=maxCategoryFdr,,drop=F]
  if(nrow(rawTable) == 0){
    return(NULL)
  }

  rawTable$URL = paste0("http://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=",rawTable$NAME)
  # For RMD
  # rawTable$NAME_URL=unlist(apply(rawTable, 1, function(x){
  #   addLinkTag(text=x['NAME'],link=x['URL'])
  # }))
  rawTable$NAME_URL=unlist(apply(rawTable, 1, function(x){
    paste0("<a href='", x['URL'], "' target='_blank'>", gsub("_", " ", x['NAME']), "</a>")
  }))
  
  descriptions=unlist(apply(rawTable, 1, function(x){
    categoryDescription=""
    cname = x['NAME']
    if(cname %in% names(gsd_map)){
      categoryDescription=gsd_map[cname]
    }else{
      gseaUrl=x['URL']
      gseaWeb<-getURL(gseaUrl)
      if(gseaWeb != ""){
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
        }else{
          return(NA)
        }
      }
    }
    return(categoryDescription)    
  }))

  rawTable$Description=descriptions

  rawTable$Core=round(rawTable$SIZE * as.numeric(str_extract(rawTable$LEADING.EDGE, "\\d+")) / 100)

  rawTable<-rawTable[,c("NAME", "NAME_URL", "Description", "SIZE", "NES", "FDR.q.val", "Core")]
  
  return(rawTable)
}

parse_gsea_subfolder=function(resultDirSub, gname, gseaCategory, is_positive){
  pstr = ifelse(is_positive, "pos", "neg")
  title = ifelse(is_positive, "Positive-regulated", "Negative-regulated")
  pattern = paste0("gsea_report_for_na_", pstr, '_\\d+\\.(tsv|xls)$')
        
  gseaTableFile<-list.files(resultDirSub,pattern=pattern,full.names=TRUE)

  if (length(gseaTableFile)!=1) {
    warning(paste("Can't find ", title, "GSEA table file in", resultDirSub))
    return(NULL)
  }
  
  rawTable<-processGseaTable(gseaTableFile)
  if(is.null(rawTable)){
    warning(paste("Can't find significant", title, "GSEA gene set in", resultDirSub))
    return(NULL)
  }
  return(rawTable)
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
  #maxCategory=0
  
  if(is_singlecell){
    gsea_files<-files
    gsea_files$Comparisons<-gsea_files$compName
    comparisons<-rev(unique(gsea_files$Comparisons))
  }else{
    gsea_files<-files
    gsea_files$Comparisons<-gsea_files$V2
    comparisons<-unique(gsea_files$Comparisons)
  }

  j=1
  for (j in 1:length(comparisons)){
    comparison<-comparisons[j]
    comp_files<-gsea_files[gsea_files$Comparisons == comparison,]

    if(print_rmd){
      cat(paste0("\n\n", gsea_prefix, "## ", comparison, "\n\n"))
    }else{
      cat(comparison, "\n")
    }
  
    has_enriched=FALSE

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

        if(!print_rmd){
          cat("  ", gseaCategory, "\n")
        }

        pos = parse_gsea_subfolder(resultDirSub, gname, gseaCategory, TRUE)
        neg = parse_gsea_subfolder(resultDirSub, gname, gseaCategory, FALSE)

        if (is.null(pos) & is.null(neg)){
          next
        }

        has_enriched = TRUE

        prefix = paste0(gname, ".", gseaCategory)

        top_enriched = NULL
        if(!is.null(pos)){
          pos_file=paste0(target_folder, prefix, ".pos.csv")
          write.csv(pos[,c(2:ncol(pos))], pos_file, row.names=F)
          result[nrow(result) + 1,] = c(comparison, gname, gseaCategory, "pos_file", pos_file)
          top_enriched = head(pos, 10)
        }

        if(!is.null(neg)){
          neg_file=paste0(target_folder, prefix, ".neg.csv")
          write.csv(neg[,c(2:ncol(neg))], neg_file, row.names=F)
          result[nrow(result) + 1,] = c(comparison, gname, gseaCategory, "neg_file", neg_file)
          top_enriched = rbind(top_enriched, head(neg, 10))
        }

        enriched_file=paste0(target_folder, prefix, ".enriched.csv")
        write.csv(top_enriched, enriched_file, row.names=F)
        result[nrow(result) + 1,] = c(comparison, gname, gseaCategory, "enriched_file", enriched_file)

        top_enriched$NES = as.numeric(top_enriched$NES)
        top_enriched$Core = as.numeric(top_enriched$Core)
        top_enriched$FDR.q.val = as.numeric(top_enriched$FDR.q.val)

        enriched_png=paste0(enriched_file, ".png")
        result[nrow(result) + 1,] = c(comparison, gname, gseaCategory, "enriched_png", enriched_png)

        top_enriched<-top_enriched[order(top_enriched$NES, decreasing=T),]
        top_enriched$NAME<-factor(top_enriched$NAME, levels=top_enriched$NAME)
        g<-ggplot(top_enriched, aes(NES, NAME)) + geom_point(aes(size=Core, color=FDR.q.val)) + geom_vline(xintercept=0) + theme_bw() + ylab("") + xlab("Normalized enrichment score")

        height=max(1200, 3000/40*nrow(top_enriched))
        ncharmax=max(nchar(as.character(top_enriched$NAME)))
        width=max(3000, ncharmax * 40 + 1000)
        ggsave(enriched_png, g, width=width, height=height, dpi=300, units="px", bg="white")

        enriched_pdf=paste0(enriched_file, ".pdf")
        ggsave(enriched_pdf, g, width=width, height=height, dpi=300, units="px", bg="white")

        if(print_rmd){
          cat(paste0("\n\n<img src='", enriched_png, "'>\n\n"))
        }
      }

      if(is_singlecell){
        break
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

save_gsea_rmd<-function(files, resFile, rmd_prefix=""){
  source("reportFunctions.R")
  if(nrow(files) == 0){
    return("# No geneset with FDR < 0.05 detected")
  }
  result<-""
  comparisons = unique(files$comparison)
  comparisons = comparisons[order(comparisons)]
  comparison<-comparisons[1]
  for(comparison in comparisons){
    result<-paste0(result, paste0("\n\n", rmd_prefix, "# ", comparison, "\n\n"))

    comp_files<-files[files$comparison==comparison,,drop=FALSE]

    categories = unique(comp_files$category)
    category=categories[1]
    for(category in categories){
      result<-paste0(result, paste0("\n\n", rmd_prefix, "## ", category, "\n\n"))

      cat_files<-comp_files[comp_files$category == category,,drop=F]
      file_map<-split(cat_files$file_path, cat_files$file_key)

      if(!is.null(file_map$pos_file)){
        result<-paste0(result, paste0("\n\n", rmd_prefix, "### Positive-regulated\n\n"))
        result<-paste0(result, getPagedTable(filepath=file_map$pos_file, row.names=0, escape=FALSE))
      }

      if(!is.null(file_map$neg_file)){
        result<-paste0(result, paste0("\n\n", rmd_prefix, "### Negative-regulated\n\n"))
        result<-paste0(result, getPagedTable(filepath=file_map$neg_file, row.names=0, escape=FALSE))
      }

      if(!is.null(file_map$enriched_png)){
        result<-paste0(result, paste0("\n\n", rmd_prefix, "### Top enriched gene sets\n\n"))
        result<-paste0(result, getFigure(file_map$enriched_png, out_width="100%"))
      }
    }
  }

  writeLines(result, resFile)
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

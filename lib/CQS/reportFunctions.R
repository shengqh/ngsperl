#duplicate with reportFunctions.Rmd

###################################
#report functions start
###################################

library(knitr)
library(reshape2)
library(ggplot2)
library(data.table)
library(DT)
library(RCurl)
library(htmltools)
library(kableExtra)
library(DT)

knitr::opts_chunk$set(echo = TRUE)

addHtmlLinkTag<-function(text,link) {
  result<-paste0("<a href='",link,"' target='_blank'>",text,"</a>")
  return(result)
}

addLinkTag<-function(text,link) {
  result<-paste0("[", text, "](", link, ")")
  return(result)
}

figRef <- local({
  tag <- numeric()
  created <- logical()
  used <- logical()
  function(label, caption, prefix = options("figcap.prefix"), 
           sep = options("figcap.sep"), prefix.highlight = options("figcap.prefix.highlight"), trunk.eval=TRUE) {
    if(trunk.eval){ 
      i <- which(names(tag) == label)
      if (length(i) == 0) {
        i <- length(tag) + 1
        tag <<- c(tag, i)
        names(tag)[length(tag)] <<- label
        used <<- c(used, FALSE)
        names(used)[length(used)] <<- label
        created <<- c(created, FALSE)
        names(created)[length(created)] <<- label
      }
      if (!missing(caption)) {
        created[label] <<- TRUE
        paste0(prefix.highlight, prefix, " ", i, sep, prefix.highlight, 
               " ", caption)
      } else {
        used[label] <<- TRUE
        paste(prefix, tag[label])
      }
    }
  }
})
options(figcap.prefix = "Figure", figcap.sep = ":", figcap.prefix.highlight = "**")


tabRef <- local({
  tag <- numeric()
  created <- logical()
  used <- logical()
  function(label, caption, prefix = options("tabcap.prefix"), 
           sep = options("tabcap.sep"), prefix.highlight = options("tabcap.prefix.highlight"), trunk.eval=TRUE) {
    if(trunk.eval){
      i <- which(names(tag) == label)
      if (length(i) == 0) {
        i <- length(tag) + 1
        tag <<- c(tag, i)
        names(tag)[length(tag)] <<- label
        used <<- c(used, FALSE)
        names(used)[length(used)] <<- label
        created <<- c(created, FALSE)
        names(created)[length(created)] <<- label
      }
      if (!missing(caption)) {
        created[label] <<- TRUE
        paste0(prefix.highlight, prefix, " ", i, sep, prefix.highlight, 
              " ", caption)
      } else {
        used[label] <<- TRUE
        paste(prefix, tag[label])
      }
    }
  }
})
options(tabcap.prefix = "Table", tabcap.sep = ":", tabcap.prefix.highlight = "**")

check_and_include_graphics<-function(graphicFile) {
  if (!is.null(graphicFile[1]) & file.exists(graphicFile[1])) {
    include_graphics(graphicFile)
  }
}

printTable<-function(filepath, row.names=1){
  if(row.names > 0){
    tbl<-data.frame(fread(filepath, check.names=F), row.names=row.names)
  }else{
    tbl<-data.frame(fread(filepath, check.names=F))
  }
  print(kable_styling(kable(tbl)))
}

printPagedTable<-function(filepath, row.names=1, escape=TRUE, digits=0, nsmall=0){
  if(row.names > 0){
    tbl<-data.frame(fread(filepath, check.names=F), row.names=row.names)
  }else{
    tbl<-data.frame(fread(filepath, check.names=F))
  }

  if(digits > 0){
    tbl <- tbl %>% dplyr::mutate_if(is.numeric, format, digits=digits, nsmall=nsmall)
  }

  DT::datatable(tbl, 
                extensions = c('FixedColumns','FixedHeader'),
                rownames = row.names > 0,
                escape = escape,
                options = list( scrollX=TRUE, 
                                paging=TRUE))
}

getPagedTable<-function(filepath, row.names=1, escape=TRUE, digits=0, nsmall=0){
  return(paste0("\n```{r,echo=FALSE,results='asis'}\nprintPagedTable('", filepath, "', ", row.names, ",", escape, ",", digits, ",", nsmall, ")\n```\n\n"))
}

getTable<-function(filepath, row.names=1){
  return(paste0("\n```{r,echo=FALSE,results='asis'}\nprintTable('", filepath, "', ", row.names, ")\n```\n\n"))
}

getFigure<-function(filepath, in_details=FALSE){
  if(in_details){
    return(paste0("\n```{r,echo=FALSE,results='asis'}\ncheck_and_include_graphics('details/", filepath, "')\n```\n\n"))
  }else{
    return(paste0("\n```{r,echo=FALSE,results='asis'}\ncheck_and_include_graphics('", filepath, "')\n```\n\n"))
  }
}

find_module_folder=function(files,pattern) {
  moduleFolders=sapply(strsplit(files[,1],"\\/"), function(x) {ind=which(x=="result");x[ind-1]})
  interestedModuleInd=grep(pattern,moduleFolders)
  return(interestedModuleInd)
}


###################################
#report functions end
###################################


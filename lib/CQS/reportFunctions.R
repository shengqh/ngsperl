#duplicate with reportFunctions.Rmd

###################################
#report functions start
###################################

copy_or_download<-function(source, target=NULL){
  if(is.null(target)){
    target = paste0("./", basename(source))
  }
  if(file.exists(source)){
    file.copy(source, target, overwrite=TRUE)
  }else{
    if(file.exists(target)){      
      file.remove(target)
    }
    download.file(source, target, "auto")
  }
}

load_install<-function(library_name, library_sources=library_name){
  if(!require(library_name, character.only = T)){
    BiocManager::install(library_sources, ask=FALSE)
  }
  library(library_name, character.only = T)
}

load_install("knitr")
load_install("reshape2")
load_install("ggplot2")
load_install("data.table")
load_install("DT")
load_install("RCurl")
load_install("htmltools")
load_install("kableExtra")
load_install("digest")
load_install("ggpubr")
load_install("patchwork")
#load_install("htmlTable")
load_install("tools")
load_install("ggExtra")
#load_install("gtsummary")
#load_install("flextable")
load_install("dplyr")
load_install("ComplexHeatmap")
load_install("grid")

options(figcap.prefix = "Figure", figcap.sep = ":", figcap.prefix.highlight = "**")
options(tabcap.prefix = "Table", tabcap.sep = ":", tabcap.prefix.highlight = "**")

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

is_file_exists<-function(filename){
  if(is.null(filename)){
    return(FALSE)
  }

  if(is.na(filename)){
    return(FALSE)
  }

  return(file.exists(filename))
}

check_and_include_graphics<-function(graphicFile) {
  if (is_file_exists(graphicFile[1])) {
    include_graphics(graphicFile)
  }
}

output_table<-function(tbl, caption=NULL, description=NULL){
  if(!is.null(caption)){
    new_caption = tabRef(caption, description)
  }else{
    new_caption = NULL
  }

  kable(tbl, caption= new_caption) %>%
    kable_styling() %>%
    htmltools::HTML()
}

print_table_from_file<-function(filepath, row.names=1, caption=NULL, description=NULL){
  if(row.names > 0){
    tbl<-data.frame(fread(filepath, check.names=F), row.names=row.names)
  }else{
    tbl<-data.frame(fread(filepath, check.names=F))
  }

  output_table(tbl, caption, description)
}

print_table<-function(tbl, round_value=3, byDT=FALSE, row.names=TRUE){
  if(round_value > 0){
    tbl <- tbl %>% dplyr::mutate(across(where(is.numeric), round, round_value))
  }
  if(byDT){
    DT::datatable(tbl, rownames = row.names, extensions = "Buttons", options = list(dom = "Bfrtip", buttons = c("excel", "csv")))
  }else{
    print(kable(tbl, row.names = row.names))
  }
}


get_table_description<-function(category, filepath, description){
  result = "\n```{r,echo=FALSE,results='asis'}\n"
  result = paste0(result, "print_table_from_file('", filepath, "', ", row.names, ",'", category, "','", description, "')\n```\n\n")
  return(result)
}

output_paged_table<-function(tbl, rownames=TRUE, escape=TRUE, digits=0, nsmall=0){
  if(digits > 0){
    tbl <- tbl %>% dplyr::mutate_if(is.numeric, format, digits=digits, nsmall=nsmall)
  }

  DT::datatable(tbl, 
                rownames = rownames, 
                escape = escape,
                extensions = "Buttons", 
                options = list( dom = "Bfrtip", 
                                buttons = c("excel", "csv")))
}

printPagedTable<-function(filepath, row.names=1, escape=TRUE, digits=0, nsmall=0){
  if(row.names > 0){
    tbl<-data.frame(fread(filepath, check.names=F), check.names=FALSE, row.names=row.names)
  }else{
    tbl<-data.frame(fread(filepath, check.names=F), check.names=FALSE)
  }

  output_paged_table( tbl=tbl,
                      rownames=row.names > 0,
                      escape=escape,
                      digits=digits, 
                      nsmall=nsmall)
}

getPagedTable<-function(filepath, row.names=1, escape=TRUE, digits=0, nsmall=0){
  return(paste0("\n```{r,echo=FALSE,results='asis'}\nprintPagedTable(filepath='", filepath, "', row.names=", row.names, ", escape=", escape, ", digits=", digits, ", nsmall=", nsmall, ")\n```\n\n"))
}

getTable<-function(filepath, row.names=1){
  return(paste0("\n```{r,echo=FALSE,results='asis'}\nprint_table_from_file('", filepath, "', ", row.names, ")\n```\n\n"))
}

getFigure<-function(filepath, in_details=FALSE, out_width=NULL){
  if(!is.null(out_width)){
    result = paste0("\n```{r,echo=FALSE,results='asis',out.width='", out_width,  "'}\n")
  }else{
    result = "\n```{r,echo=FALSE,results='asis'}\n"
  }
  if(in_details){
    return(paste0(result, "check_and_include_graphics('details/", filepath, "')\n```\n\n"))
  }else{
    return(paste0(result, "check_and_include_graphics('", filepath, "')\n```\n\n"))
  }
}

getFigure_width_height<-function(filepath, in_details=FALSE, fig.width=NULL, fig.height=NULL){
  if(!is.null(fig.width)){
    if(!is.null(fig.height)){
      result = paste0("\n```{r,echo=FALSE,results='asis',fig.width=", fig.width, ",fig.height=", fig.height, "}\n")
    }else{
      result = paste0("\n```{r,echo=FALSE,results='asis',fig.width=", fig.width, "}\n")
    }
  }else{
    if(!is.null(fig.height)){
      result = paste0("\n```{r,echo=FALSE,results='asis',fig.height=", fig.height, "}\n")
    }else{
      result = paste0("\n```{r,echo=FALSE,results='asis'}\n")
    }
  }
  if(in_details){
    return(paste0(result, "check_and_include_graphics('details/", filepath, "')\n```\n\n"))
  }else{
    return(paste0(result, "check_and_include_graphics('", filepath, "')\n```\n\n"))
  }
}

get_figure_description<-function(category, filepath, description){
  return(paste0("```{r,echo=FALSE,results='asis', fig.align='center', fig.cap=figRef('", category, "', '",gsub("_", " ", description), "', trunk.eval=file.exists('", filepath, "'))}\n",
"  check_and_include_graphics('", filepath, "')\n```\n"))
}

find_module_folder=function(files,pattern) {
  moduleFolders=sapply(strsplit(files[,1],"\\/"), function(x) {ind=which(x=="result");x[ind-1]})
  interestedModuleInd=grep(pattern,moduleFolders)
  return(interestedModuleInd)
}

get_date_str = function(){
  format(Sys.time(), "%Y%m%d")
}

check_md5<-function(filepath, expect_md5, return_md5=FALSE){
  if(!file.exists(filepath)){
    stop("File not exists: ", filepath)
  }
  md5=tools::md5sum(filepath)

  if(expect_md5 == ""){
    if(return_md5){
      return(md5)
    }else{
      cat(basename(filepath), "md5=", md5, "\n")
    }
  }else{
    if(md5 != expect_md5){
      stop("md5 not match, expect ", expect_md5, " but got ", md5, " for file ", filepath)
    }
  }
}

theme_rotate_x_axis_label <- function() {
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

theme_bw3 <- function (axis.x.rotate=F) { 
  result = theme_bw() +
    theme(
      strip.background = element_rect(fill = NA, colour = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),			
      axis.line = element_line(colour = "black", linewidth = 0.5),
      plot.title = element_text(hjust = 0.5)
    )
  if (axis.x.rotate){
    result = result + theme_rotate_x_axis_label()
  }
  
  return(result)
}

get_hist_density<-function(data, x, title=x, bins=20){
  ggplot(data, aes(x=!!sym(x))) + geom_histogram(aes(y = ..density..), colour = 1, fill = "white", bins=bins) + geom_density() + ggtitle(title) + theme_bw3()
}

show_descriptive_statistics<-function(data){
  dd = data
  dd$fakevar = 1
  dd$fakevar <- factor(dd$fakevar, levels = c(1), labels = c("Subject"))
  #label(dd$fakevar) <- "Subject"

  dd_formula = paste0(paste0(colnames(data), collapse=" + "), " ~ fakevar")
  print_descriptive_statistics(as.formula(dd_formula), dd, test = FALSE)
}

print_descriptive_statistics<-function(formula, data, test = TRUE, overall = FALSE, continuous = 5, ...){
  library(Hmisc)
  output <- summaryM(formula = formula,
                      data = data, 
                      test = test, 
                      overall = overall, 
                      continuous = continuous, ...)
  latex_tbl = latex(output, html=TRUE, width=0.8 )
  cat(latex_tbl)
}

factor_by_count<-function(vec){
  tbl=table(vec)
  tbl=tbl[tbl > 0]
  tbl=tbl[order(tbl, decreasing=T)]
  res=factor(vec, levels=names(tbl))
  return(res)
}

get_log_cpm<-function(counts, prefix=NULL, filterCPM=TRUE, transform=TRUE){
  library(edgeR)
  
  dge <- DGEList(counts)
  dge <- calcNormFactors(dge)

  cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = F)
  if(!is.null(prefix)){
    write.csv(cpm, paste0(prefix, ".cpm.csv"))
  }

  if(filterCPM){
    ncpm1=rowSums(cpm >= 1)
    keep=ncpm1 > floor(ncol(cpm)/2)
    dge <- dge[keep, ,keep.lib.sizes=TRUE]
  }

  logcpm <- cpm(dge, normalized.lib.sizes = TRUE, log = T, prior.count = 2)
  if(!is.null(prefix)){
    write.csv(logcpm, paste0(prefix,".logcpm.csv"))
  }

  if(transform){
    logcpmt <- t(logcpm)
    return(logcpmt)
  }else{
    return(logcpm)
  }
}

read_file_map<-function(file_list_path, sep="\t", header=F, do_unlist=TRUE){
  tbl=fread(file_list_path, header=header, data.table=FALSE)

  result<-split(tbl$V1, tbl$V2)
  if(do_unlist){
    result<-unlist(result)
  }
  return(result)
}

get_option_value<-function(myoptions, key, required=FALSE, default=NULL){
  if(key %in% names(myoptions)){
    return(myoptions[[key]])
  }
  
  if(required){
    stop("Missing required option: ", key)
  }

  return(default)
}

get_legend_width <- function(g, by="max") {
  gg <- ggplotGrob(g)
  
  # Find the legend element and extract its width
  legend_grob <- gg$grobs[[which(sapply(gg$grobs, function(x) {x$name}) == "guide-box")]]
  if(by=="max"){
    legend_width <- ceiling(convertWidth(max(legend_grob$width), "in", valueOnly = TRUE))
  }else{
    legend_width <- ceiling(convertWidth(sum(legend_grob$width), "in", valueOnly = TRUE))
  }

  return(legend_width)
}

get_freq_table<-function(df, column){
df |> 
  dplyr::count(!!sym(column)) |> 
  dplyr::arrange(desc(n)) |>
  dplyr::rename("count"="n")
}

cluster_dotplot<-function(gdata, column1="features.plot", column2="id", value_column="avg.exp.scaled", dim="both"){
  load_install("textshape")

  mdata=acast(gdata, as.formula(paste0(column1, "~", column2)), value.var=value_column)

  mdata=cluster_matrix(mdata, dim=dim, method="ward.D2")
  gdata[,column1]=factor(gdata[,column1], levels=rownames(mdata))
  gdata[,column2]=factor(gdata[,column2], levels=colnames(mdata))

  return(gdata)
}

summary_tableby = function(dat, formula, test=F, total=T) {
  library(arsenal)
  mycontrols  <- tableby.control( test=test, 
                                  total=total,
                                  numeric.stats=c("Nmiss","median","q1q3","range"),
                                  cat.stats=c("Nmiss","countpct"),
                                  stats.labels=list(Nmiss='Missing', median='Median', range='Range'),
                                  digits=1, 
                                  digits.p=2, 
                                  digits.pct=0)

  result = arsenal::tableby(formula, data = dat, control=mycontrols)
  summary(result, title = "", width = 3, pfootnote = TRUE)  
}

###################################
#report functions end
###################################


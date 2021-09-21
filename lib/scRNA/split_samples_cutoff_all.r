
source("split_samples_utils.r")

library(Seurat)
library(ggplot2)

#devtools::install_github("shengqh/cutoff")
#install.packages("bbmle")
library(zoo)
library(reshape2)
library(gridExtra)
library(ggExtra)

library(bbmle)
library(choisycutoff)

#still keep it here in case we need it later.
my_startval <- function(values,D1="normal",D2="normal",cutoff_point=0) {
  if(cutoff_point >0){
    thresh = cutoff_point
  }else{
    den <- tryCatch(
      expr = {
        density(values, bw="SJ")
      },
      error = function(e){ 
        density(values)
      }
    )
    w=1
    x=den$x
    y=den$y
    y.smooth=den$y
    n <- length(y.smooth)
    y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
    delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
    i.max <- which(delta <= 0) + w
    res=data.frame(x=x[i.max], i=i.max, y=y[i.max])
    res=res[res$x>0,]
    res=res[order(res$y, decreasing = T),]

    if(nrow(res)>2){
      res=res[1:2,]
    }

    xx=x[x>min(res$x) & x<max(res$x)]
    yy=y[x>min(res$x) & x<max(res$x)]
    yy.min=min(yy)
    ii.min=which(yy==yy.min)
    thresh=xx[ii.min]
  }

  sel <- values<thresh
  data1 <- values[sel]
  data2 <- values[!sel]
  lambda <- length(data1)/length(values)
  param1 <- MASS::fitdistr(data1,D1)$est
  param2 <- MASS::fitdistr(data2,D2)$est
  out <- c(param1,param2,lambda)
  names(out) <- c("mu1","sigma1","mu2","sigma2","lambda")
  return(out)
}

get_cutoff<-function(values, prefix, cutoff_startval=NA){
  my_out <- em(values,"normal","normal", cutoff_startval=cutoff_startval)
  cut_off <- cutoff(my_out)
  png(paste0(prefix, ".cutoff.png"), width=2000, height=1600, res=300)
  hist(values,200,F,xlab="concentration",ylab="density", main=NULL,col="grey")
  lines(density(values),lwd=1.5,col="blue")
  lines(my_out,lwd=1.5,col="red")
  abline(v=cut_off[1],lwd=1.5,col="brown")
  abline(v=my_out$start.val,lwd=1.5,col="green")
  dev.off()
  return(list(em_out=my_out, cut_off=cut_off))
}

files_lines=read.table(parSampleFile1, sep="\t")
files=split(files_lines$V1, files_lines$V2)

if(parSampleFile2 != ''){
  cutoffs_lines=read.table(parSampleFile2, sep="\t")
  cutoffs=split(cutoffs_lines$V1, cutoffs_lines$V2)
}else{
  cutoffs=list()
}

params_lines=read.table(parSampleFile3, sep="\t")
params=split(params_lines$V1, params_lines$V2)
params$hto_ignore_exists=ifelse(params$hto_ignore_exists=="0", FALSE, TRUE)

idx=14
for(idx in c(1:length(files))){
  fname=names(files)[idx]
  output_prefix = paste0(fname, ".HTO")
  output_file=paste0(output_prefix, ".csv")
  
  #if(file.exists(output_file) & params$hto_ignore_exists){
  #  next
  #}

  rdsfile=files[[idx]]
  cat(fname, ":", rdsfile, " ...\n")

  obj=read_hto(rdsfile, output_prefix)
  
  cutoff_point=ifelse(fname %in% names(cutoffs), as.numeric(cutoffs[[fname]]), NA)

  tagnames=rownames(obj[["HTO"]])
  data <- FetchData(object=obj, vars=tagnames)
  write.csv(data, file=paste0(output_prefix, ".data.csv"))
  
  res=NULL
  tagname=tagnames[1]  
  for (tagname in tagnames) {
    values=data[,tagname]
    values=values[values>0]
    cat(paste0("get cutoff of ", tagname, " ...\n"))
    my_cutoff=get_cutoff(values, paste0(output_prefix, "_", tagname), cutoff_point)
    data[,paste0(tagname,"_pos")] = ifelse(data[,tagname]>my_cutoff$cut_off[1], tagname, "Negative")
    res<-rbind(res, data.frame(
      start.val=my_cutoff$em_out$start.val,
      mu1=my_cutoff$em_out$param["mu1"], 
      sigma1=my_cutoff$em_out$param["sigma1"], 
      mu2=my_cutoff$em_out$param["mu2"], 
      sigma2=my_cutoff$em_out$param["sigma2"], 
      cutoff=my_cutoff$cut_off[1]))
  }
  rownames(res)=tagnames
  write.csv(res, paste0(output_prefix, ".estimated.csv"))
  
  cat(paste0("get classification ...\n"))
  class_names=paste0(tagnames, "_pos")
  data$HTO_classification=unlist(apply(data, 1, function(x){
    xx=unique(x[class_names])
    if (length(xx) > 1){
      xx = xx[xx != "Negative"]
      if (length(xx) > 1) {
        return("Doublet")
      }
    }
    return(xx)
  }))
  
  data$HTO_classification.global=unlist(apply(data, 1, function(x){
    xx=unique(x[class_names])
    if (length(xx) > 1){
      xx = xx[xx != "Negative"]
      if (length(xx) > 1) {
        return("Doublet")
      }else{
        return("Singlet")
      }
    }
    return("Negative")
  }))
  
  obj[["HTO_classification"]] = data$HTO_classification
  obj[["HTO_classification.global"]] = data$HTO_classification.global
  
  output_post_classification(obj, output_prefix)
}

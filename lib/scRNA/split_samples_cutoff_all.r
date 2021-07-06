source("split_samples_utils.r")
library(Seurat)
library(ggplot2)

#devtools::install_github("choisy/cutoff")
#install.packages("bbmle")
library(choisycutoff)
library(zoo)
library(reshape2)
library(gridExtra)
library(ggExtra)

my_startval <- function(values,D1="normal",D2="normal",cutoff_point=0) {
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
  #the highest peaks should be the negative one, positive one is always in the right side.
  res=res[res$x>=res$x[1],]

  if(cutoff_point > 0){
    if(res$x[1] > cutoff_point){
      res=rbind(res[1,,drop=F], res[res$x < cutoff_point,,drop=F])
    }else{
      res=rbind(res[1,,drop=F], res[res$x > cutoff_point,,drop=F])
    }
  }
  
  if(nrow(res)>2){
    if(nrow(res)>2){
      res=res[1:2,]
    }
  }

  xx=x[x>min(res$x) & x<max(res$x)]
  yy=y[x>min(res$x) & x<max(res$x)]
  yy.min=min(yy)
  ii.min=which(yy==yy.min)
  thresh=xx[ii.min]
  
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

t=1e-64
my_em<-function(values, data_name="em", D1="normal", D2="normal", t=1e-64, cutoff_point=0){
  start <- as.list(my_startval(values, D1, D2, cutoff_point))
  
  D1b <- choisycutoff:::hash[[D1]]
  D2b <- choisycutoff:::hash[[D2]]
  lambda0 <- 0
  with(start, {
    while (abs(lambda0 - mean(lambda)) > t) {
      lambda <- mean(lambda)
      lambda0 <- lambda
      distr1 <- lambda * D1b(values, mu1, sigma1)
      distr2 <- (1 - lambda) * D2b(values, mu2, sigma2)
      lambda <- distr1/(distr1 + distr2)
      mLL2 <- function(mu1, sigma1, mu2, sigma2) return(choisycutoff:::mLL(mu1, 
                                                                     sigma1, mu2, sigma2, lambda, values, D1b, D2b))
      start <- as.list(log(c(mu1 = mu1, sigma1 = sigma1, 
                             mu2 = mu2, sigma2 = sigma2)))
      out <- bbmle::mle2(mLL2, start, "Nelder-Mead")
      coef <- out@coef
      coef_n <- names(coef)
      names(coef) <- NULL
      for (i in 1:4) assign(coef_n[i], exp(coef[i]))
    }
    out <- list(lambda = lambda, param = exp(out@coef), D1 = D1, 
                D2 = D2, deviance = out@min, data = values, data_name = data_name, 
                out = out, t = t)
    class(out) <- "em"
    return(out)
  })
}

my_cutoff<-function (object, t = 1e-64, nb = 10, distr = 2, type1 = 0.05, level = 0.95) 
{
  coef <- object$out@coef
  the_names <- names(coef)
  coef <- exp(mc2d::rmultinormal(nb, coef, as.vector(object$out@vcov)))
  coef <- as.list(data.frame(t(coef)))
  coef <- lapply(coef, function(x) {
      names(x) <- the_names
      return(as.list(x))
  })
  out <- sapply(coef, function(x) choisycutoff:::lci0(x, mean(object$lambda), 
      choisycutoff:::hash[[object$D1]], choisycutoff:::hash[[object$D2]], object$data, t))
  lambda <- rnorm(nb, out[1, ], out[2, ])
  coef <- sapply(coef, function(x) unlist(x))
  the_names <- c(rownames(coef), "lambda")
  coef <- rbind(coef, lambda)
  coef <- lapply(as.data.frame(coef), function(x) {
      names(x) <- the_names
      return(as.list(x))
  })
  out <- sapply(coef, function(x) with(x, choisycutoff:::cutoff0(mu1, 
      sigma1, mu2, sigma2, lambda, object$D1, object$D2, distr, type1)))
  out <- MASS::fitdistr(out, "normal")
  the_mean <- out$estimate["mean"]
  level <- (1 - level)/2
  level <- c(level, 1 - level)
  ci <- the_mean + qt(level, Inf) * out$sd["mean"]
  out <- c(the_mean, ci)
  names(out) <- c("Estimate", paste(100 * level, "%"))
  return(out)
}

get_cutoff<-function(values, prefix, cutoff_point=0){
  my_out <- my_em(values,"normal","normal", cutoff_point=cutoff_point)
  cut_off <- my_cutoff(my_out)
  png(paste0(prefix, ".cutoff.png"), width=2000, height=1600, res=300)
  hist(values,200,F,xlab="concentration",ylab="density", main=NULL,col="grey")
  lines(density(values),lwd=1.5,col="blue")
  lines(my_out,lwd=1.5,col="red")
  abline(v=cut_off[1],lwd=1.5,col="brown")
  dev.off()
  return(cut_off[1])
}

files_lines=read.table(parSampleFile1, sep="\t")
files=split(files_lines$V1, files_lines$V2)

cutoffs_lines=read.table(parSampleFile2, sep="\t")
cutoffs=split(cutoffs_lines$V1, cutoffs_lines$V2)

params_lines=read.table(parSampleFile3, sep="\t")
params=split(params_lines$V1, params_lines$V2)
params$hto_ignore_exists=ifelse(params$hto_ignore_exists=="0", FALSE, TRUE)

allres=NULL
idx=3
for(idx in c(1:length(files))){
  fname=names(files)[idx]
  output_prefix = paste0(fname, ".HTO")
  output_file=paste0(output_prefix, ".csv")
  allres<-rbind(allres, data.frame(File=output_file, Sample=fname))
  
  if(file.exists(output_file) & params$hto_ignore_exists){
    next
  }

  h5file=files[[idx]]
  cat(fname, ":", h5file, " ...\n")

  obj=read_hto(h5file, output_prefix, params$hto_regex)
  
  cutoff_point=ifelse(fname %in% names(cutoffs), as.numeric(cutoffs[[fname]]), 0)

  tagnames=rownames(obj[["HTO"]])
  data <- FetchData(object=obj, vars=tagnames)
  
  tagname=tagnames[1]  
  for (tagname in tagnames) {
    values=data[,tagname]
    values=values[values>0]
    cat(paste0("get cutoff of ", tagname, " ...\n"))
    cutoff=get_cutoff(values, paste0(output_prefix, "_", tagname), cutoff_point)
    data[,paste0(tagname,"_pos")] = ifelse(data[,tagname]>cutoff, tagname, "Negative")
  }
  
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

write.csv(allres, file=paste0(outFile, ".htofile.csv"), row.names=F)

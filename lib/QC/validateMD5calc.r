rm(list=ls()) 
outFile='veteran.csv'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/ravi_shah_projects/shengq2/20230327_rnaseq_veteran_hg38/md5_validation/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
options = read_file_map('fileList3.txt', do_unlist=FALSE)
threads = to_numeric(options$threads, 1)

files = read.table("fileList1.txt", header=F, sep="\t")
files$name = basename(files$V1)
colnames(files) = c("path", "sample", "name")

md5_files = read.table("fileList2.txt", header=F, sep="\t")

md5_values_list = lapply(md5_files$V1, function(x) { 
  return(read.table(x, header=F, sep=" ")) 
})

md5_values = do.call(rbind, md5_values_list)[,c(1,3)]
colnames(md5_values) = c("md5_expect", "name")

df = merge(x = files, y = md5_values, by = "name", all.x = TRUE)

do_md5<-function(path){
  cat("md5sum", path, "...\n")
  tools::md5sum(path)
}
md5_calc = parallel::mclapply(df$path, do_md5, mc.cores=threads)

df$md5_calc = unlist(md5_calc)
df$md5_check = df$md5_calc == df$md5_expect

write.csv(df, file=outFile, row.names=F)
failed = df[!df$md5_check,]

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')

if(nrow(failed) > 0){
  write.csv(failed, file="failed.csv", row.names=F)
  stop("md5 check failed, please check failed.csv for details")
}else{
  if(file.exists("failed.csv")){
    file.remove("failed.csv")
  }
  cat("md5 check passed\n")
}

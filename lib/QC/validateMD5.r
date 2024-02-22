rm(list=ls()) 
outFile='veteran.csv'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/nobackup/h_cqs/ravi_shah_projects/shengq2/20230327_rnaseq_veteran_hg38/md5_merge/result/veteran.md5.txt'
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/ravi_shah_projects/shengq2/20230327_rnaseq_veteran_hg38/md5_validation/result')

### Parameter setting end ###

files = read.table(parSampleFile1, header=F, sep="\t")
files$name = unlist(lapply(files$V1, basename))
colnames(files) = c("path", "sample", "name")

md5_files = read.table("fileList2.txt", header=F, sep="\t")
md5_values_list = lapply(md5_files$V1, function(x) { 
  return(read.table(x, header=F, sep=" ")) 
})
md5_values = do.call(rbind, md5_values_list)[,c(1,3)]
colnames(md5_values) = c("md5_expect", "name")

md5_calc = read.table(parFile1, sep=" ", header=F)
md5_calc = md5_calc[,c(1,3)]
names(md5_calc) = c("md5_calc", "name")
md5_calc$name=basename(md5_calc$name)

df = merge(x = files, y = md5_calc, by = "name", all.x = TRUE)
df = merge(x = df, y = md5_values, by = "name", all.x = TRUE)
df$md5_check = ifelse(df$md5_calc == df$md5_expect, "OK", "FAIL")

write.csv(df, file=outFile, row.names=F)

failed = subset(df, md5_check == "FAIL")

if(nrow(failed) > 0){
  write.csv(failed, file="failed.csv", row.names=F)
  stop("md5 check failed, please check failed.csv for details")
}else{
  if(file.exists("failed.csv")){
    file.remove("failed.csv")
  }
  cat("md5 check passed\n")
}

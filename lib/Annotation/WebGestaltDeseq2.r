library("rmarkdown")

annoFiles<-read.table(parSampleFile1, header=F, sep="\t", stringsAsFactors = F)
deseq2Files<-read.table(parSampleFile2, header=F, sep="\t", stringsAsFactors = F)

comparisons<-unique(annoFiles$V2)
comparison<-comparisons[1]
for (comparison in comparisons){
  compAnnoFiles<-annoFiles[annoFiles$V2 == comparison,]$V1
  compDeseq2File<-deseq2Files[deseq2Files$V2 == comparison,1]
  
  compAnnoFile<-compAnnoFiles[1]
  for(compAnnoFile in compAnnoFiles){
    plotData <- list(WebGestaltFile = compAnnoFile,
                     deseq2File = compDeseq2File)
    
    output_path <- paste0(normalizePath(compAnnoFile), ".html")
    
    output_dir = "E:/temp"
    output_file = "temp1.html"
    #output_dir = dirname(output_path)
    #output_file = basename(output_path)
    
    cat("Output report to:", output_path, "\n")
    rmarkdown::render("WebGestaltDeseq2.rmd",
                      output_dir = output_dir,
                      output_file = output_file,
                      params = list(data = plotData))
  }
}

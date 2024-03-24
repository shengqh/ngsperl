rm(list=ls()) 
outFile='iSGS_cell_atlas'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20240320_scRNA_iSGS_cell_atlas/cellbender_expect_cell/result')

### Parameter setting end ###

source("scRNA_func.r")
source("reportFunctions.R")

get_num_cells<-function(fileList){
  flist<-fread(fileList, header=FALSE)
  ncells<-unlist(lapply(flist$V1, function(x){
    cat("Reading ", x, " ...\n")
    obj<-read_scrna_data(x)
    return(ncol(obj$counts))
  }))
  return(data.frame(sample=flist$V2, num_cells=ncells))
}

tbl=get_num_cells(parSampleFile1)
write.csv(tbl, paste0(outFile, ".num_cells.csv"), row.names=FALSE)

tbl$perl<-paste0("    '", tbl$sample, "' => ['", tbl$num_cells, "'],")
perlines=c("  cellbender_expect_count => {", tbl$perl, "  },")
writeLines(perlines, con=paste0(outFile, ".cellranger_count_perl.txt"))

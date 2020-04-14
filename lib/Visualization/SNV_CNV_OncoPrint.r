
options(bitmapType='cairo')

library(ComplexHeatmap)

inputFiles = parSampleFile1

options <- read.table(parSampleFile2, sep="\t", stringsAsFactors = F, header=F)
rownames(options)<-options$V2

width = as.numeric(options["picture_width", "V1"])
height = as.numeric(options["picture_height", "V1"])
sampleNamePattern=options["sampleNamePattern", "V1"]

MISSENSE_color=options["MISSENSE_color", "V1"]
MISSENSE_height=as.numeric(options["MISSENSE_height", "V1"])
TRUNC_color=options["TRUNC_color", "V1"]
TRUNC_height=as.numeric(options["TRUNC_height", "V1"])
DUP_color=options["DUP_color", "V1"]
DUP_height=as.numeric(options["DUP_height", "V1"])
DEL_color=options["DEL_color", "V1"]
DEL_height=as.numeric(options["DEL_height", "V1"])

cnvFile = parFile1
cnvdata<-read.delim(cnvFile, as.is=T, header=T, sep="\t", stringsAsFactors = F)

inputFileData<-read.delim(inputFiles, header=F, sep="\t", stringsAsFactors = F)

inputFile<-inputFileData$V1[[1]]
for(inputFile in inputFileData$V1){
  oncoprint<-read.delim(inputFile,as.is=T,header=TRUE,sep="\t", row.names=1, stringsAsFactors = F)
  oncoprint[oncoprint==" "]<-NA
  oncoprint[oncoprint==""]<-NA
  
  for(idx in c(1:nrow(cnvdata))){
    sample = cnvdata[idx, "File"]
    gene = cnvdata[idx, "Feature"]
    cnv = cnvdata[idx, "CNV"]
    oldvalue = oncoprint[gene, sample]
    if(is.na(oldvalue)){
      oncoprint[gene, sample] = paste0(cnv, ";")
    }else{
      oncoprint[gene, sample] = paste0(oldvalue, cnv, ";")
    }
  }


  outputTextFile = paste0(outputDirectory, "/", basename(inputFile), ".snv_cnv.txt")
  write.csv(oncoprint, file=outputTextFile, quote=F, na="")
  
  col = c("MISSENSE" = MISSENSE_color, "TRUNC" = TRUNC_color, "DUP" = DUP_color, "DEL" = DEL_color)
  
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    MISSENSE = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*MISSENSE_height, gp = gpar(fill = col["MISSENSE"], col = NA))
    },
    TRUNC = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*TRUNC_height, gp = gpar(fill = col["TRUNC"], col = NA))
    },
    DUP = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*DUP_height, gp = gpar(fill = col["DUP"], col = NA))
    },
    DEL = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*DEL_height, gp = gpar(fill = col["DEL"], col = NA))
    }
  )
  
  if(width == 0){
    width=max(2000, ncol(oncoprint) * 70 + 300)
  }

  if(height == 0){
    height=max(1000, nrow(oncoprint) * 70 + 300)
  }

  ##oncoprint
  png(paste0(outputTextFile, ".png"), width=width, height=height, res=300)
  ht=oncoPrint(oncoprint, get_type = function(x) strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, 
               column_title = "",
               show_column_names = T,
               right_annotation = rowAnnotation(
                 rbar = anno_oncoprint_barplot(
                   width = unit(1, "cm"))),
               heatmap_legend_param = list(title = "Genetic alternations", at = c("MISSENSE", "TRUNC", "DUP", "DEL"), 
                                           labels = c("Missense mutation", "Truncating mutation     ", "CNV duplication", "CNV deletion")))
  draw(ht)
  dev.off()
}


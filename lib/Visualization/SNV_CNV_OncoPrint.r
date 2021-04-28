
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

MANUAL_order=options["MANUAL_order", "V1"] != "0"

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

cnvFile = parFile1
cnvdata<-read.delim(cnvFile, as.is=T, header=T, sep="\t", stringsAsFactors = F,check.names=F)

inputFileData<-read.delim(inputFiles, header=F, sep="\t", stringsAsFactors = F)

inputFile<-inputFileData$V1[[1]]
for(inputFile in inputFileData$V1){
  oncoData<-read.delim(inputFile,as.is=T,header=TRUE,sep="\t", row.names=1, stringsAsFactors = F,check.names=F)
  oncoData[oncoData==" "]<-NA
  oncoData[oncoData==""]<-NA
  
  for(idx in c(1:nrow(cnvdata))){
    sample = cnvdata[idx, "File"]
    gene = cnvdata[idx, "Feature"]
    cnv = cnvdata[idx, "CNV"]
    oldvalue = oncoData[gene, sample]
    if(is.na(oldvalue)){
      oncoData[gene, sample] = paste0(cnv, ";")
    }else{
      oncoData[gene, sample] = paste0(oldvalue, cnv, ";")
    }
  }
  oncoData[is.na(oncoData)]<-''
  variantCount = apply(oncoData !='', 1, sum)
  variantOrder = order(variantCount, decreasing = T)
  
  df<-data.frame("V1"=colnames(oncoData), "V2"="")
  if (parSampleFile3 != ""){
    sg<-read.table(parSampleFile3, sep="\t", header=F, stringsAsFactors=F)
    df<-rbind(df, sg)
  }
  
  for (name in unique(df$V2)){
    cellSamples = df$V1[df$V2 == name]
    curOncoData = oncoData[, cellSamples,drop=F]
    
    if (name == ""){
      outputTextFile = paste0(outputDirectory, "/", basename(inputFile), ".snv_cnv.txt")
    }else{
      outputTextFile = paste0(outputDirectory, "/", basename(inputFile), ".", name, ".snv_cnv.txt")
    }
    write.csv(curOncoData, file=outputTextFile, quote=F, na="")
    
    if(width == 0){
      width=max(2000, ncol(curOncoData) * 70 + 300)
    }
    
    if(height == 0){
      height=max(1000, nrow(curOncoData) * 70 + 300)
    }
    
    if (MANUAL_order){
      columnOrder = c(1:length(cellSamples))
    } else{
      columnOrder = NULL
    }
    
    ##oncoprint
    ht=oncoPrint(curOncoData, get_type = function(x) strsplit(x, ";")[[1]],
                 alter_fun = alter_fun, col = col, 
                 column_title = "",
                 column_order = columnOrder,
                 row_order = variantOrder,
                 show_column_names = T,
                 right_annotation = rowAnnotation(
                   rbar = anno_oncoprint_barplot(
                     width = unit(1, "cm"))),
                 heatmap_legend_param = list(title = "Genetic alternations", at = c("MISSENSE", "TRUNC", "DUP", "DEL"), 
                                             labels = c("Missense mutation", "Truncating mutation     ", "CNV duplication", "CNV deletion")))
    png(paste0(outputTextFile, ".png"), width=width, height=height, res=300)
    draw(ht)
    dev.off()
    pdf(paste0(outputTextFile, ".pdf"), width=10, height=10*height/width)
    draw(ht)
    dev.off()
  }
}


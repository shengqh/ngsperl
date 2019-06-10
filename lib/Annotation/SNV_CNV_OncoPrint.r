#rm(list=ls()) 
#outFile='linton_exomeseq_3321'
#parSampleFile1=''
#parSampleFile2=''
#parSampleFile3=''
#parSampleFile4=''
#parFile1='/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/bwa_refine_gatk4_hc_gvcf_vqsr_filterMAF_annovar_filter/result/linton_exomeseq_3321.freq0.001.snv.missense.oncoprint.tsv'
#parFile2='/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/GATK4_CNV_Germline_11_AnnotationGenesPlot/result/linton_exomeseq_3321.position.txt.slim'
#parFile3=''


options(bitmapType='cairo')

library(ComplexHeatmap)

inputFiles = parSampleFile1
cnvFile = parFile1

cnvdata<-read.table(cnvFile, as.is=T, header=T, sep="\t")

inputFileData<-read.table(inputFiles, sep="\t", stringsAsFactors = F)

inputFile<-inputFileData$V1[[1]]
for(inputFile in inputFileData$V1){
  oncoprint<-read.table(inputFile,as.is=T,header=TRUE,sep="\t", row.names=1)
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
  
  col = c("MISSENSE" = "darkgreen", "TRUNC" = "red", "DUP" = "brown", "DEL" = "blue")
  
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    MISSENSE = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.80, gp = gpar(fill = col["MISSENSE"], col = NA))
    },
    TRUNC = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.60, gp = gpar(fill = col["TRUNC"], col = NA))
    },
    DUP = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.40, gp = gpar(fill = col["DUP"], col = NA))
    },
    DEL = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.20, gp = gpar(fill = col["DEL"], col = NA))
    }
  )
  
  width=max(2000, ncol(oncoprint) * 70 + 300)
  height=1600
  
  ##oncoprint
  outputFile = paste0(outputDirectory, "/", basename(inputFile), ".snv_cnv.png")
  png(outputFile, width=width, height=height, res=300)
  ht=oncoPrint(oncoprint, get_type = function(x) strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, 
               column_title = "",
               show_column_names = T,
               row_barplot_width = unit(0.5, "cm"),
               heatmap_legend_param = list(title = "Genetic alternations", at = c("MISSENSE", "TRUNC", "DUP", "DEL"), 
                                           labels = c("Missense mutation", "Truncating mutation", "CNV Dup", "CNV Del")))
  draw(ht)
  dev.off()
}


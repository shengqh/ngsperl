options(bitmapType='cairo')

library(ComplexHeatmap)

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0){
  setwd('/scratch/cqs/fanr1/202111_ICARE_BEST/DNA/20211207_DNAseq_Analysis/bwa_g4_refine_gatk4_SNV_05_filter/result')
  inputFile = "human_exomeseq.freq0.001.snv.missense.tsv"
  outputFile = "human_exomeseq.freq0.001.snv.missense.top10.oncoprint.tsv"
  optionFile = "onco_options.txt"
  #genelist = c("BRAF","RAS", "NTRK2", "CDKN2A", "CDKN2B", "NF1", "KMT2D", "RB1", "MMR", "ARID2", "ATM")
  genelist = NA
}else{
  inputFile = args[1]
  outputFile = args[2]
  optionFile = args[3]
  genelist = ifelse(length(args) >= 4, args[4:length(args)], NA)
}

cat("inputFile=", inputFile, "\n")
cat("outputFile=", outputFile, "\n")
cat("optionFile=", optionFile, "\n")

options <- read.table(optionFile, sep="\t", stringsAsFactors = F, header=F)
rownames(options)<-options$V2

width = as.numeric(options["picture_width", "V1"])
height = as.numeric(options["picture_height", "V1"])
sampleNamePattern=options["sampleNamePattern", "V1"]

BACKGROUND_color=options["BACKGROUND_color", "V1"]
MISSENSE_color=options["MISSENSE_color", "V1"]
MISSENSE_height=as.numeric(options["MISSENSE_height", "V1"])
TRUNC_color=options["TRUNC_color", "V1"]
TRUNC_height=as.numeric(options["TRUNC_height", "V1"])

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = BACKGROUND_color, col = NA))
  },
  MISSENSE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.66, gp = gpar(fill = MISSENSE_color, col = NA))
  },
  TRUNC = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = TRUNC_color, col = NA))
  }
)
col = c("MISSENSE" = MISSENSE_color, "TRUNC" = TRUNC_color)

mutdata<-read.delim(inputFile,as.is=T,header=TRUE,sep="\t",stringsAsFactors = F)
mutdata[is.na(mutdata)] <- 0
cnames=colnames(mutdata)

countNotNumIndex<-which(grepl("Format", cnames))
if (length(countNotNumIndex)==0) {
  index<-1;
  indecies<-c()
} else {
  index<-max(countNotNumIndex)+1
  indecies<-c(1:(index-1))
}
cnames<-cnames[index:length(cnames)]

samples<-cnames[grep(sampleNamePattern,cnames)]
if(length(samples) == 0){
  stop(paste0("No sample matches the pattern ", sampleNamePattern))
}

if(is.na(genelist)){
  mdata<-mutdata[,c('Gene.refGene', cnames)]
  library(dplyr)
  mdata2 <- aggregate(. ~ Gene.refGene, data = mdata, max)
  rownames(mdata2)<-mdata2$Gene.refGene
  mdata2<-mdata2[,c(2:ncol(mdata2)),drop=F]
  mdata3<-ifelse(mdata2 > 0, 1, 0)
  mdata4<-apply(mdata3, 1, sum)
  mdata5<-mdata4[order(mdata4,decreasing=T)]
  if(length(mdata5) >= 10){
    genelist<-names(mdata5)[1:10]
  }else{
    genelist<-names(mdata5)
  }
  cat("top genes: ", genelist, "\n")
}

geneind<-mutdata$Gene.refGene %in% genelist

mutdata_gene_samples<-mutdata[geneind,c("Func.refGene", "Gene.refGene", "ExonicFunc.refGene", samples)]
sample_start<-4
sample_end<-ncol(mutdata_gene_samples)

##oncoprint
oncoprint<-NULL
gene=genelist[3]
for (gene in genelist){
  cat("gene=", gene, "\n")
  genedata=mutdata_gene_samples[mutdata_gene_samples$Gene.refGene == gene,]
  
  value = c()
  i=69
  for (i in sample_start:sample_end){
    sample<-colnames(mutdata_gene_samples)[i]
    mut_sample<-genedata[,i]
    if (sum(mut_sample)>0) {
      curvalues=c()
      for (j in 1: length(mut_sample)) {
        if (genedata[j,i]>0) {
          type<-genedata[j,3]
          if (type=="frameshift deletion" | type=="frameshift insertion" | type=="frameshift substitution" | type=="stopgain" | type=="stoploss" ) 
            type="TRUNC"
          else if (type=="nonsynonymous SNV") 
            type="MISSENSE"
          else if (genedata[j,1]=="splicing") 
            type="TRUNC"
          else
            #stop(paste0("unknown genotype ", type))
            type="UNKNOWN"
          curvalues=c(curvalues, paste0(type, ";"))
        } 
      }
      curvalue=paste0(unique(curvalues), collapse = '')
    }else {
      curvalue=" "
    }
    value=c(value, curvalue)
  }
  if(any(value != " ")){
    oncoprint<-rbind(oncoprint,c(gene,value))
  }
}
colnames(oncoprint)<-c("Gene",colnames(mutdata_gene_samples)[sample_start:sample_end])
write.table(oncoprint,file=outputFile,quote=F,row.names=F,sep="\t")

rownames(oncoprint)=oncoprint[,1]
oncoprint=oncoprint[,c(2:ncol(oncoprint)),drop=F]

if(width == 0){
  width=max(2000, ncol(oncoprint) * 70 + 300)
}

if(height == 0){
  height=max(1000, nrow(oncoprint) * 70 + 300)
}

png(paste0(outputFile, ".png"), width=width, height=height, res=300)
ht=oncoPrint(oncoprint, get_type = function(x) strsplit(x, ";")[[1]],
             alter_fun = alter_fun, col = col, 
             column_title = "",
             show_column_names = T,
             right_annotation = rowAnnotation(
               rbar = anno_oncoprint_barplot(
                 width = unit(1, "cm"))),
             heatmap_legend_param = list(title = "Mutation", at = c("MISSENSE", "TRUNC"), 
                                         labels = c("Missense  ", "Truncating  ")))
draw(ht)
dev.off()

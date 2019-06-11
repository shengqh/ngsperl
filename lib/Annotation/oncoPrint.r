options(bitmapType='cairo')

library(ComplexHeatmap)

args = commandArgs(trailingOnly=TRUE)
inputFile = args[1]
outputFile = args[2]
width = args[3]
height = args[4]
sampleNamePattern=args[5]
genelist=args[6:length(args)]

DEBUG = FALSE
if (DEBUG) {
  inputFile = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/bwa_refine_gatk4_hc_gvcf_vqsr_filterMAF_annovar_filter/result/linton_exomeseq_3321.freq0.001.snv.missense.tsv"
  outputFile = "/scratch/cqs/shengq2/macrae_linton/20190517_linton_exomeseq_3321_human/bwa_refine_gatk4_hc_gvcf_vqsr_filterMAF_annovar_filter/result/linton_exomeseq_3321.freq0.001.snv.missense.oncoprint.tsv"
  width = 0
  height = 0
  sampleNamePattern="."
  genelist=c("LDLR")
}

cat("inputFile=", inputFile, "\n")
cat("outputFile=", outputFile, "\n")
cat("width=", width, "\n")
cat("height=", height, "\n")
cat("sampleNamePattern=", sampleNamePattern, "\n")
cat("genelist=", genelist, "\n")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  MISSENSE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.66, gp = gpar(fill = "darkgreen", col = NA))
  },
  TRUNC = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "red", col = NA))
  }
)
col = c("MISSENSE" = "darkgreen", "TRUNC" = "red")

mutdata<-read.delim(inputFile,as.is=T,header=TRUE,sep="\t",stringsAsFactors = F)
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

geneind<-mutdata$Gene.refGene %in% genelist

mutdata_gene_samples<-mutdata[geneind,c("Func.refGene", "Gene.refGene", "ExonicFunc.refGene", samples)]
sample_start<-4
sample_end<-ncol(mutdata_gene_samples)

##oncoprint
oncoprint<-NULL
gene=genelist[length(genelist)]
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
          if (type=="frameshift deletion" | type=="frameshift insertion" | type=="frameshift substitution" | type=="stopgain") 
            type="TRUNC"
          else if (type=="nonsynonymous SNV") 
            type="MISSENSE"
          else if (genedata[j,1]=="splicing") 
            type="TRUNC"
          else
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
oncoprint=oncoprint[,c(2:ncol(oncoprint))]

if(width == 0){
  width=max(2000, ncol(oncoprint) * 70 + 300)
}

if(height == 0){
  height=max(1600, nrow(oncoprint) * 70 + 300)
}

png(paste0(outputFile, ".png"), width=width, height=height, res=300)
ht=oncoPrint(oncoprint, get_type = function(x) strsplit(x, ";")[[1]],
             alter_fun = alter_fun, col = col, 
             column_title = "",
             show_column_names = T,
             row_barplot_width = unit(0.5, "cm"),
             heatmap_legend_param = list(title = "Genetic alternations", at = c("MISSENSE", "TRUNC"), 
                                         labels = c("Missense mutation", "Truncating mutation")))
draw(ht)
dev.off()

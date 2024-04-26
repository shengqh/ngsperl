options(bitmapType='cairo')

library(ComplexHeatmap)

args = commandArgs(trailingOnly=TRUE)
inputFile = args[1]
outputPrefix = args[2]
sampleNamePattern=args[3]
genelist=args[4:length(args)]

if(!all(is.na(genelist))){
  if(any(grepl(",", genelist))){
    genelist = gsub(" ", "", genelist)
    genelist = sort(unlist(strsplit(genelist, ",")))
  }
}

#setwd("/scratch/cqs/shengq1/dnaseq/20161013_liuqi_gene_panel/bwa_refine_hc_gvcf_vqsr_annovar_filter/result/")
#inputFile = "Adenoma.exac0.001.tsv"
#outputPrefix = "Adenoma.exac0.001"
#sampleNamePattern="cur_"
#genelist=c("APC", "ROBO2", "KRAS", "NRAS", "BRAF", "TP53")

cat("inputFile=", inputFile, "\n")
cat("outputPrefix=", outputPrefix, "\n")
cat("sampleNamePattern=", sampleNamePattern, "\n")
cat("genelist=", genelist, "\n")

oncoFile = paste0(outputPrefix, ".oncoprinter.txt")
mutationmapperFile = paste0(outputPrefix, ".mutationmapper.txt")

mutdata<-read.table(inputFile,as.is=T,header=TRUE,sep="\t")
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

mutdata_gene_samples<-mutdata[geneind,c("Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene", samples)]
sample_start<-5
sample_end<-ncol(mutdata_gene_samples)

##oncoprint
oncoprint<-NULL
i=69
for (i in sample_start:sample_end){
	sample<-colnames(mutdata_gene_samples)[i]
	mut_sample<-mutdata_gene_samples[,i]
	if (sum(mut_sample)>0) {
		for (j in 1: length(mut_sample)) {
			if (mut_sample[j]>0) {
				gene=mutdata_gene_samples$Gene.refGene[j]
				
				info<-mutdata_gene_samples$AAChange.refGene[j]
				info<-unlist(strsplit(info,":"))
				pchange<-info[length(info)]
				
				type<-mutdata_gene_samples$ExonicFunc.refGene[j]
				if (type=="frameshift deletion" | type=="frameshift insertion" | type=="frameshift substitution" | type=="stopgain"){ 
					type="TRUNC"
				}else if (type=="nonsynonymous SNV"){ 
					type="MISSENSE"
				}else if (mutdata_gene_samples$Func.refGene[j]=="splicing"){ 
					type="TRUNC"
					pchage="splice"
				}else{
					type="UNKNOWN"
				}
				
				oncoprint<-rbind(oncoprint,c(sample,gene, pchange,type))
			} 
		}
	}
}
colnames(oncoprint)=c("Sample","Gene","Alteration", "Type")
write.table(oncoprint,file=oncoFile,quote=F,row.names=F,sep="\t")

###mutationmapper
mutationmapper<-NULL
for (i in sample_start:sample_end){
	sample<-colnames(mutdata_gene_samples)[i]
	mut_sample<-mutdata_gene_samples[,i]
	if (sum(mut_sample)>0) {
		for (j in 1: length(mut_sample)) {
			if (mut_sample[j]>0) {
				if (mutdata_gene_samples$Func.refGene[j]!="splicing"){
					gene=mutdata_gene_samples$Gene.refGene[j]
					
					info<-mutdata_gene_samples$AAChange.refGene[j]
					info<-unlist(strsplit(info,":"))
					pchange<-info[length(info)]
					
					type<-mutdata_gene_samples$ExonicFunc.refGene[j]
					if (type=="frameshift deletion" | type=="frameshift substitution") type="Frame_Shift_Del"
					else if (type=="frameshift insertion") type="Frame_Shift_Ins"
					else if (type=="stopgain") type="Nonsense_Mutation"
					else if (type=="nonsynonymous SNV") type="Missense_Mutation"
					
					mutationmapper<-rbind(mutationmapper,c(gene, sample,pchange, type))
				} 
			} 
		}
	}
}
colnames(mutationmapper)<-c("Hugo_Symbol","Sample_ID","Protein_Change","Mutation_Type")
write.table(mutationmapper,file=mutationmapperFile,row.names=F,quote=F,sep="\t")

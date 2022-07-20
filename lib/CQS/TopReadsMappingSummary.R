rm(list=ls()) 
outFile='RA_4893_2'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/scratch/vickers_lab/projects/20220707_4893_2_RA_smRNA_mouse_v5_forSpcount/class_independent/identical_sequence_count_table/result/RA_4893_2_sequence.read.count'
parFile2=''
parFile3=''
readFilesModule=c('Host miRNA','Host tRNA','Host snRNA','Host snoRNA','Host rRNA','Host other small RNA','Host Genome','Refseq Bacteria','Microbiome Bacteria','Environment Bacteria','Fungus','Algae','Virus','Non host tRNA','Non host rRNA'); textSize=9;groupTextSize=10;

setwd('/gpfs23/scratch/vickers_lab/projects/20220707_4893_2_RA_smRNA_mouse_v5_forSpcount/data_visualization/top_sequence_mapped_in_categories/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
options(bitmapType='cairo')

resultFile<-outFile
top100ReadFile<-parFile1
readFileList<-parSampleFile1
groupFileList<-parSampleFile2
groupVisLayoutFileList<-parSampleFile3

top100Table<-read.delim(top100ReadFile,header=T,row.names=1,as.is=T,check.names=F)
readFiles<-read.delim(readFileList,header=F,as.is=T)[,1]
mappingTable<-NULL
for (readFile in readFiles) {
	readsTable<-read.delim(readFile,header=T,row.names=1,as.is=T,check.names=F)
	mappingResult<-rep("N",nrow(top100Table))
	names(mappingResult)<-row.names(top100Table)
	
	temp<-intersect(row.names(top100Table),row.names(readsTable))
	mappingResult[temp]<-"Y"
	mappingTable<-cbind(mappingTable,mappingResult)
}
colnames(mappingTable)<-readFilesModule
temp<-apply(mappingTable,1,function(x) length(which(x=="Y")))
mappingTable<-cbind(mappingTable,NumberOfModulesMapped=temp)

result<-cbind(mappingTable,top100Table)
write.csv(result,paste0(resultFile,".ReadsMapping.Summary.csv"))

library(VennDiagram)
library(tools)

draw_venn<-function(df, cat_sets, filename){
  x = list()
  for(cname in cat_sets){
    sequences = rownames(df)[df[,cname] == 'Y']
    x[[cname]] = sequences
  }
  venn.diagram(x, filename, imagetype=file_ext(filename))
}

if("Refseq Bacteria" %in% colnames(result)){
	draw_venn(df=result, cat_sets=c("Refseq Bacteria", "Microbiome Bacteria", "Environment Bacteria"), filename=paste0(resultFile, ".bacteria.venn.png"))

	draw_venn(df=result, cat_sets=c("Refseq Bacteria", "Fungus", "Algae", "Virus"), filename=paste0(resultFile, "nonhost_genome.venn.png"))

	draw_venn(df=result, cat_sets=c("Refseq Bacteria", "Non host tRNA", "Non host rRNA"), filename=paste0(resultFile, "nonhost_genome_library.venn.png"))

	logfiles = list.files(pattern = paste0(resultFile, ".*.log"))
	for (logfile in logfiles){
		unlink(logfile)
	}
}

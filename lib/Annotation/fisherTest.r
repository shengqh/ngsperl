options(bitmapType='cairo')

args = commandArgs(trailingOnly=TRUE)

inputFile = args[1]
outputFile = args[2]
fisherConditionFile=args[3]

#inputFile = "/scratch/cqs/shengq1/dnaseq/20161013_liuqi_gene_panel/bwa_refine_hc_gvcf_vqsr_annovar_filter_base1/result/Adenoma.base1.exac0.01.snv.tsv"
#outputFile = "/scratch/cqs/shengq1/dnaseq/20161013_liuqi_gene_panel/bwa_refine_hc_gvcf_vqsr_annovar_filter_base1/result/Adenoma.base1.exac0.01.snv.fisher.tsv"
#fisherConditionFile = "/scratch/cqs/shengq1/dnaseq/20161013_liuqi_gene_panel/documents/Base1RecurNorecur.txt"

cat("inputFile=", inputFile, "\n")
cat("outputFile=", outputFile, "\n")
cat("fisherConditionFile=", fisherConditionFile, "\n")

data=read.table(inputFile, sep="\t", header=T)

if("Chr" == colnames(data)[1]){
	rownames(data)=paste(data$Chr, data$Start, data$End, data$Ref, data$Alt, sep="_")
	infocolumns=c(1,2,3,4,5,7,9)
}else{
	rownames(data)=data$Gene
	infocolumns=c(1,2)
}

condition=read.table(fisherConditionFile, sep="\t", header=T, stringsAsFactor=F)
condition=condition[condition[,1] %in% colnames(data),]

uniqConds = unlist(unique(condition[,2]))
if(length(uniqConds) != 2){
  stop(paste0("Only two condition allowed for fisher test, now is " , paste(uniqConds, collapse=",")))
}

conditionData=data[,condition[,1]]

cond1=condition[,2]==uniqConds[1]
cond2=!cond1
cond1Count=sum(cond1)
cond2Count=sum(cond2)

ft <- apply(conditionData, 1, function(x){
			cond1data=x[cond1]
			cond2data=x[cond2]
			cond1true=sum(cond1data)
			cond1false=cond1Count-cond1true
			cond2true=sum(cond2data)
			cond2false=cond2Count-cond2true
			m = matrix(c(cond1true, cond1false, cond2true, cond2false),
					nrow=2,
					dimnames=list(cond1=c("yes", "no"),
					              cond2=c("yes", "no")))
			cat(m, "\n")
			ft=fisher.test(m)
			return(c(cond1true, cond1false, cond2true, cond2false, ft$p.value ))
		})

ft<-t(ft)
colnames(ft)<-c(paste0(uniqConds[1], c("True", "False")), paste0(uniqConds[2], c("True", "False")), "pValue" )
ft<-data.frame(ft)
ft<-ft[order(ft$pValue),]
fdata<-data.frame(data[rownames(ft), infocolumns], ft)
write.table(fdata, file=outputFile, sep="\t", row.names=F, quote=F)

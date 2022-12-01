options(bitmapType='cairo')

args = commandArgs(trailingOnly=TRUE)
positionFile = args[1]
file=tools::file_path_sans_ext(positionFile)
print(file)

library(reshape2)
library(ggplot2)
library(dplyr)

snRNAGrouping<-function(x) {
	if (all(grepl("^RNU|^RNVU|^U\\d+",head(as.character(x))))) {
		return(1) #snRNA
	} else if (all(grepl("tRNA-",head(as.character(x))))) {
		return(2) #tRNA
	} else {
		return(0)
	}
}

snRnaName2Group<-function(x,groupSnRNA=1) {
	if (groupSnRNA==1) { #snRNA
		snRnaGroup<-sapply(strsplit(as.character(x),"-|:"),function(y) y[1])
		snRnaGroup<-gsub("RNVU","RNU",snRnaGroup)
		snRnaGroup<-gsub("([A-Z]+[0-9]+)[A-Z]+","\\1",snRnaGroup)
	} else if (groupSnRNA==2) { #tRNA
		snRnaGroup<-sapply(strsplit(as.character(x),"-"),function(y) y[(grep("tRNA",y)[1]+1)])
	} else {
		snRnaGroup<-x
	}
	return(snRnaGroup)
}


fp=read.table(positionFile, sep="\t", header=T)

featureCount<- as.data.frame(fp %>% group_by(Feature) %>% summarise(max = max(Count)))
rownames(featureCount)<-featureCount$Feature
featureCount$max<-round(featureCount$max, 0)
fp$MaxCount<-featureCount[fp$Feature,"max"]
fp$Feature<-paste0(fp$Feature,"(",fp$MaxCount,"):",fp$Strand)

features=as.character(unique(fp$Feature))
if(length(features) > 100){
  features = features[1:100]
  fp=fp[fp$Feature %in% features,]
}
fp$Feature<-factor(fp$Feature, levels=rev(unique(fp$Feature)))

groupSnRNA<-snRNAGrouping(features)

fp$AbsCount<-fp$Count*fp$Percentage
fp$Percentage[which(fp$Percentage==0)]<-NA
fp$snRNAGroup<-snRnaName2Group(fp$Feature,groupSnRNA)

height=max(length(unique(fp$Feature))*100,3000)
width=height
png(paste0(file,".png"), width=width, height=height, res=300)
p<-ggplot(fp,aes(x=Position,y=Feature,size=Percentage,colour=Percentage))+
		geom_point()+
		scale_size_continuous(range = c(0.1,2))+
		scale_colour_gradient(low="indianred1",high="darkred")+
		xlim(c(-30, 120))+ 
		theme_bw() +
		theme(legend.position="none")+
		theme(text = element_text(size=25))
if (groupSnRNA) {
	p<-p+facet_grid(snRNAGroup~.,space = "free",scale="free")
}
print(p)
dev.off()

if (groupSnRNA) {
	counts=lapply(features, function(x){
				idx=which(fp$Feature==x)
				unlist(fp$Count[idx])[1]
			})
	counts=round(unlist(counts))
	names(counts)<-features
	snRNAGroup<-snRnaName2Group(features,groupSnRNA)
	names(snRNAGroup)<-features
	snRNAGroupTotalCounts<-tapply(counts,snRNAGroup[names(counts)],sum)
	fpBySnRNAGroup<-aggregate(x = fp, by = list(fp$snRNAGroup, fp$Position), FUN = function(x) if(is.numeric(x)| is.integer(x)) {sum(x)} else {x[1]})
	fpBySnRNAGroup$GroupPercentage<-fpBySnRNAGroup$AbsCount/snRNAGroupTotalCounts[fpBySnRNAGroup$snRNAGroup]
	fpBySnRNAGroup$Position<-fpBySnRNAGroup$Group.2
	
	height=max(length(unique(fpBySnRNAGroup$snRNAGroup))*100,1500)
	width=height
	
	if(groupSnRNA==1){
	  name = ".snRNAGroup"
	}else{
	  name = ".tRNAGroup"
	}
	
	png(paste0(file, name, ".png"), width=width, height=height, res=300)
	p<-ggplot(fpBySnRNAGroup,aes(x=Position,y=snRNAGroup,size=GroupPercentage,colour=GroupPercentage))+
			geom_point()+
			scale_size_continuous(range = c(0.1,2))+
			scale_colour_gradient(low="indianred1",high="darkred")+
			xlim(c(-30, 120))+ 
			theme_bw() +
			theme(legend.position="none")
	print(p)
	dev.off()
}





#############################
#Vis for all Non Host Reads: Group 1, 2, 4; tRNA; two rRNA Categry;
#############################

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")
resultFile<-outFile
readFileList<-parSampleFile1
groupFileList<-parSampleFile2
groupVisLayoutFileList<-parSampleFile3
totalCountFile<-parFile3

categoriesNames<-c("Microbiome","Environment","Fungus","tRNA","rRNAL","rRNAS")

readsMappingNames<-list()
readsMappingTable<-NULL
readFiles<-read.delim(readFileList,header=F,as.is=T)
for (i in 1:nrow(readFiles)) {
	temp<-read.delim(readFiles[i,1],header=T,row.names=1,as.is=T)
	readsMappingNames[[i]]<-row.names(temp)
	readsMappingTable<-rbind(readsMappingTable,temp)
}

######################################
#Reads Overlap: Reads were found in how many categories? For reads only in one category, which one?
######################################
readsMappingNamesTable<-table(unlist(readsMappingNames))
dataForPlot<-NULL
#maxReadCategory<-max(readsMappingNamesTable)
#for (i in 2:maxReadCategory) {
#	temp<-colSums(readsMappingTable[names(readsMappingNamesTable)[which(readsMappingNamesTable==i)],])
#	dataForPlot<-rbind(dataForPlot,temp)
#}
#row.names(dataForPlot)<-paste0("Reads Mapped to ",2:maxReadCategory," Categories")

readsInOneCategory<-names(readsMappingNamesTable)[which(readsMappingNamesTable==1)]
for (i in 1:length(readsMappingNames)) {
	temp<-intersect(readsInOneCategory,readsMappingNames[[i]])
	temp<-colSums(readsMappingTable[temp,])
	dataForPlot<-rbind(dataForPlot,temp)
}
temp<-colSums(readsMappingTable[names(readsMappingNamesTable)[which(readsMappingNamesTable>1)],])
dataForPlot<-rbind(dataForPlot,temp)
#row.names(dataForPlot)[(nrow(dataForPlot)-length(categoriesNames)+1):nrow(dataForPlot)]<-categoriesNames
row.names(dataForPlot)<-c(paste0(categoriesNames," Only"),"Mapped to more than one Category")


#Pie chart for all samples
ggpieToFile(dataForPlot,fileName=paste0(resultFile,".Piechart.png"),maxCategory=maxCategory,textSize=textSize)

#Barplot for all samples
tableBarplotToFile(dataForPlot,fileName=paste0(resultFile,".Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart
ggpieGroupToFile(dataForPlot,fileName=paste0(resultFile,".Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)


######################################
#Reads Overlap: Venn for 5 categories
######################################
library("VennDiagram")
venn.diagram1<-function (x, count=NULL,filename, height = 3000, width = 3000, resolution = 500, 
		units = "px", compression = "lzw", na = "stop", main = NULL, 
		sub = NULL, main.pos = c(0.5, 1.05), main.fontface = "plain", 
		main.fontfamily = "serif", main.col = "black", main.cex = 1, 
		main.just = c(0.5, 1), sub.pos = c(0.5, 1.05), sub.fontface = "plain", 
		sub.fontfamily = "serif", sub.col = "black", sub.cex = 1, 
		sub.just = c(0.5, 1), category.names = names(x), force.unique = TRUE,
		fill=NA,
		...) 
{
	if (is.null(count)) {
		countFun<-function(x) length(x)
	} else {
		countFun<-function(x) sum(count[x])
	}
	if (is.na(fill)) {
		if (length(x)==5) {
			fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
		} else if (length(x)==4) {
			fill = c("dodgerblue", "goldenrod1",  "seagreen3", "orchid3")
		} else if (length(x)==3) {
			fill = c("dodgerblue", "goldenrod1", "seagreen3")
		} else if (length(x)==2) {
			fill = c("dodgerblue", "goldenrod1")
		}
	}
	if (force.unique) {
		for (i in 1:length(x)) {
			x[[i]] <- unique(x[[i]])
		}
	}
	if ("none" == na) {
		x <- x
	}
	else if ("stop" == na) {
		for (i in 1:length(x)) {
			if (any(is.na(x[[i]]))) {
				stop("NAs in dataset", call. = FALSE)
			}
		}
	}
	else if ("remove" == na) {
		for (i in 1:length(x)) {
			x[[i]] <- x[[i]][!is.na(x[[i]])]
		}
	}
	else {
		stop("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"")
	}
	if (0 == length(x) | length(x) > 5) {
		stop("Incorrect number of elements.", call. = FALSE)
	}
	if (1 == length(x)) {
		list.names <- category.names
		if (is.null(list.names)) {
			list.names <- ""
		}
		grob.list <- VennDiagram::draw.single.venn(area = length(x[[1]]), 
				category = list.names, ind = FALSE,fill=fill, ...)
	}
	else if (2 == length(x)) {
		grob.list <- VennDiagram::draw.pairwise.venn(area1 = length(x[[1]]), 
				area2 = length(x[[2]]), cross.area = length(intersect(x[[1]], 
								x[[2]])), category = category.names, ind = FALSE, 
				fill=fill,
				...)
	}
	else if (3 == length(x)) {
		A <- x[[1]]
		B <- x[[2]]
		C <- x[[3]]
		list.names <- category.names
		nab <- intersect(A, B)
		nbc <- intersect(B, C)
		nac <- intersect(A, C)
		nabc <- intersect(nab, C)
		grob.list <- VennDiagram::draw.triple.venn(area1 = length(A), 
				area2 = length(B), area3 = length(C), n12 = length(nab), 
				n23 = length(nbc), n13 = length(nac), n123 = length(nabc), 
				category = list.names, ind = FALSE, list.order = 1:3, 
				fill=fill,
				...)
	}
	else if (4 == length(x)) {
		A <- x[[1]]
		B <- x[[2]]
		C <- x[[3]]
		D <- x[[4]]
		list.names <- category.names
		n12 <- intersect(A, B)
		n13 <- intersect(A, C)
		n14 <- intersect(A, D)
		n23 <- intersect(B, C)
		n24 <- intersect(B, D)
		n34 <- intersect(C, D)
		n123 <- intersect(n12, C)
		n124 <- intersect(n12, D)
		n134 <- intersect(n13, D)
		n234 <- intersect(n23, D)
		n1234 <- intersect(n123, D)
		grob.list <- VennDiagram::draw.quad.venn(area1 = length(A), 
				area2 = length(B), area3 = length(C), area4 = length(D), 
				n12 = length(n12), n13 = length(n13), n14 = length(n14), 
				n23 = length(n23), n24 = length(n24), n34 = length(n34), 
				n123 = length(n123), n124 = length(n124), n134 = length(n134), 
				n234 = length(n234), n1234 = length(n1234), category = list.names, 
				ind = FALSE, fill=fill,...)
	}
	else if (5 == length(x)) {
		A <- x[[1]]
		B <- x[[2]]
		C <- x[[3]]
		D <- x[[4]]
		E <- x[[5]]
		list.names <- category.names
		n12 <- intersect(A, B)
		n13 <- intersect(A, C)
		n14 <- intersect(A, D)
		n15 <- intersect(A, E)
		n23 <- intersect(B, C)
		n24 <- intersect(B, D)
		n25 <- intersect(B, E)
		n34 <- intersect(C, D)
		n35 <- intersect(C, E)
		n45 <- intersect(D, E)
		n123 <- intersect(n12, C)
		n124 <- intersect(n12, D)
		n125 <- intersect(n12, E)
		n134 <- intersect(n13, D)
		n135 <- intersect(n13, E)
		n145 <- intersect(n14, E)
		n234 <- intersect(n23, D)
		n235 <- intersect(n23, E)
		n245 <- intersect(n24, E)
		n345 <- intersect(n34, E)
		n1234 <- intersect(n123, D)
		n1235 <- intersect(n123, E)
		n1245 <- intersect(n124, E)
		n1345 <- intersect(n134, E)
		n2345 <- intersect(n234, E)
		n12345 <- intersect(n1234, E)
		grob.list <- VennDiagram::draw.quintuple.venn(area1 = countFun(A), 
				area2 = countFun(B), area3 = countFun(C), area4 = countFun(D), 
				area5 = countFun(E), n12 = countFun(n12), n13 = countFun(n13), 
				n14 = countFun(n14), n15 = countFun(n15), n23 = countFun(n23), 
				n24 = countFun(n24), n25 = countFun(n25), n34 = countFun(n34), 
				n35 = countFun(n35), n45 = countFun(n45), n123 = countFun(n123), 
				n124 = countFun(n124), n125 = countFun(n125), n134 = countFun(n134), 
				n135 = countFun(n135), n145 = countFun(n145), n234 = countFun(n234), 
				n235 = countFun(n235), n245 = countFun(n245), n345 = countFun(n345), 
				n1234 = countFun(n1234), n1235 = countFun(n1235), n1245 = countFun(n1245), 
				n1345 = countFun(n1345), n2345 = countFun(n2345), n12345 = countFun(n12345), 
				category = list.names, ind = FALSE,fill=fill, ...)
	}
	else {
		stop("Invalid size of input object")
	}
	if (!is.null(sub)) {
		grob.list <- add.title(gList = grob.list, x = sub, pos = sub.pos, 
				fontface = sub.fontface, fontfamily = sub.fontfamily, 
				col = sub.col, cex = sub.cex)
	}
	if (!is.null(main)) {
		grob.list <- add.title(gList = grob.list, x = main, pos = main.pos, 
				fontface = main.fontface, fontfamily = main.fontfamily, 
				col = main.col, cex = main.cex)
	}
	grid.newpage()
	grid.draw(grob.list)
	return(1)
#	return(grob.list)
}
if (length(readsMappingNames)>5) {
	dataForPlot<-readsMappingNames[1:5]
} else {
	dataForPlot<-readsMappingNames
}
names(dataForPlot)<-categoriesNames[1:5]
for (i in 1:ncol(readsMappingTable)) {
	reads2count<-readsMappingTable[,i]
	names(reads2count)<-row.names(readsMappingTable)
	png(paste0(resultFile,".",colnames(readsMappingTable)[i],".venn.png"),res=300,height=2000,width=2000)
	venn.diagram1(dataForPlot,count=reads2count)
	dev.off()
}



###################################################
#Group1 and Group2 reads mapping table overlap
###################################################
temp1<-intersect(readsMappingNames[[1]],readsMappingNames[[2]])
temp2<-setdiff(readsMappingNames[[1]],temp1)
temp3<-setdiff(readsMappingNames[[2]],temp1)

temp1<-colSums(readsMappingTable[temp1,])
temp2<-colSums(readsMappingTable[temp2,])
temp3<-colSums(readsMappingTable[temp3,])
dataForPlot<-rbind(BothCategories=temp1,MicrobiomeOnly=temp2,EnvironmentOnly=temp3)

#Pie chart for all samples
ggpieToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Piechart.png"),maxCategory=maxCategory,textSize=textSize)

#Barplot for all samples
tableBarplotToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Barplot.png"),totalCountFile=totalCountFile,maxCategory=maxCategory,textSize=textSize)

#Group Pie chart
ggpieGroupToFile(dataForPlot,fileName=paste0(resultFile,".MicrobiomeVsEnvironment.Group.Piechart.png"),groupFileList=groupFileList,
		outFileName=paste0(resultFile,".PercentGroups.csv"),maxCategory=maxCategory,textSize=groupTextSize,visLayoutFileList=groupVisLayoutFileList)



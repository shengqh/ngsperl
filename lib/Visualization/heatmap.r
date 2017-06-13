options(bitmapType='cairo')

args <- commandArgs(trailingOnly=TRUE)

#args<-c(
#"/scratch/cqs/shengq1/brown/20170410_heatmap/result/TNF/TNF.gff.filelist",
#"blue",
#"TNF",
#"TRUE",

#getting mapped tables
gffFiles<-read.table(args[1], header=F, stringsAsFactor=F, sep="\t")

#set the plot color
plot_color= args[2]

#set the plot name
plot_name= args[3]

#set rpm
rpm_string = args[4]
if(rpm_string == "TRUE"){
	rpm=TRUE
}else{
	rpm= FALSE
}
print('RPM is:')
print(rpm_string)

plotHeatmaps <- function(mapped_data,plot_color,plot_path,names_list,rpm=TRUE){
	if(rpm == TRUE){
		y_title = 'rpm/bp'
	}else{
	  y_title = 'raw reads/bp'
	}
			
	#setting up color scaling
	colorSpectrum <- colorRampPalette(c("white",plot_color))(100)

  minq <- function(x) { quantile(x, ra.rm=TRUE,prob=0.6) }
	minValue = min(unlist(lapply(mapped_data, minq)))
	print(paste('min value is', minValue))

  maxq <- function(x) { quantile(x, ra.rm=TRUE,prob=0.95) }
	maxValue = max(unlist(lapply(mapped_data, maxq)))
	print(paste('max value is', maxValue))
	
	color_cuts <- seq(minValue,maxValue,length=100)
	
	trueMin = min(unlist(lapply(mapped_data, min)))
	trueMax = max(unlist(lapply(mapped_data, max)))
	
	color_cuts <- c(trueMin, color_cuts,trueMax)
	colorSpectrum <- c(colorSpectrum[1],colorSpectrum)

	referenceOrder = order(apply(mapped_data[[1]],1,mean,na.rm=TRUE))	
	
	png(filename=plot_path,width = 1500,height = 1600)
	layout(matrix(data=c(1:length(mapped_data)),nrow=1))
	for(mapIndex in c(1:length(mapped_data))){
    data_matrix=mapped_data[[mapIndex]]
    image(1:ncol(data_matrix),1:nrow(data_matrix),t(data_matrix[referenceOrder,]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",xlab="",ylab="",main=names_list[mapIndex])
  }
	dev.off()	
}

mapped_tables=lapply(gffFiles$V1, read.delim)
trans <- function(x) { as.matrix(x[,3:ncol(x)]) }
mapped_data<-lapply(mapped_tables, trans)

plot_path = paste0(plot_name,'.unscaled_heat.png')
print(plot_path)
plotHeatmaps(mapped_data,plot_color,plot_path,gffFiles$V2,rpm)

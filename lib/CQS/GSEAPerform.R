# TODO: Add comment
# 
# Author: zhaos
###############################################################################



runGSEA<-function(preRankedGeneFile,resultDir=NULL,gseaJar="/home/zhaos/bin/gsea-3.0.jar",gseaDb="/scratch/TBI/Data/Reference_genome/GSEA/human/V6.0",
#		gseaCategories=c("h.all.v6.0.symbols.gmt","c2.all.v6.0.symbols.gmt","c5.all.v6.0.symbols.gmt","c6.all.v6.0.symbols.gmt","c7.all.v6.0.symbols.gmt"),
		gseaCategories=c("h.all.v6.0.symbols.gmt"),
		gseaReportTemplt="GSEAReport.Rmd",
		makeReport=FALSE
)
{
	fileToName=c("h"="HallmarkGeneSets","c1"="PositionalGeneSets","c2"="CuratedGeneSets","c3"="MotifGeneSets","c4"="ComputationalGeneSets","c5"="GOGeneSets","c6"="OncogenicGeneSets","c7"="ImmunologicGeneSets")
	
	if (is.null(resultDir)) {
		gesaResultDir<-paste0(preRankedGeneFile,".gsea")
	} else {
		gesaResultDir<-paste0(resultDir,"/",basename(preRankedGeneFile),".gsea")
	}
	if (file.exists(gesaResultDir)) {
		warning(paste0(gesaResultDir," folder exists! Will delete all files in it and regenerate GSEA results."))
		unlink(gesaResultDir, recursive = TRUE)
	}
	
	for (gseaCategory in gseaCategories) {
		gseaCategoryName=strsplit(gseaCategory,"\\.")[[1]][1]
		if (gseaCategoryName %in% names(fileToName)) {
			gseaCategoryName<-fileToName[gseaCategoryName]
		}
		runCommond=paste0("java -Xmx8198m -cp ",gseaJar," xtools.gsea.GseaPreranked -gmx ",gseaDb,"/",gseaCategory,
				"-collapse false -nperm 1000 -rnk ",preRankedGeneFile," -scoring_scheme weighted -make_sets true -rpt_label '",gseaCategoryName,"' -plot_top_x 20 -set_max 500 -set_min 15 -out ",gesaResultDir)
		print(runCommond)
		system(runCommond)
	}
	
	if (makeReport) {
		library('rmarkdown')
		rmarkdown::render(gseaReportTemplt,output_file=paste0(basename(preRankedGeneFile),".gsea.html"),output_dir=resultDir)
		
	}
}


preRankedGeneFileTable=parSampleFile1
preRankedGeneFileTable=read.delim(preRankedGeneFileTable,header=F,as.is=T)

if(!exists("gseaCategories")){
  gseaCategories<-"h.all.v6.0.symbols.gmt"
}

for (i in 1:nrow(preRankedGeneFileTable)) {
	preRankedGeneFile=preRankedGeneFileTable[i,1]
	compName=preRankedGeneFileTable[i,2]
	
	resultDir=getwd()
	runGSEA(preRankedGeneFile,resultDir=resultDir,makeReport=TRUE,gseaJar=gseaJar,gseaDb=gseaDb,gseaCategories=gseaCategories)
	
}





###############################################################################
# Author: Shilin Zhao, Quanhu Sheng
###############################################################################

runGSEA<-function(preRankedGeneFile,resultDir=NULL,gseaJar="gsea-cli.sh",gseaDb="/scratch/cqs_share/references/gsea/v7.0",
#		gseaCategories=c('h.all.v7.0.symbols.gmt', 'c2.all.v7.0.symbols.gmt', 'c5.all.v7.0.symbols.gmt', 'c6.all.v7.0.symbols.gmt', 'c7.all.v7.0.symbols.gmt'),
		gseaCategories=c("h.all.v7.0.symbols.gmt"),
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

		if (grepl("cli.sh", gseaJar)){
			runCommand=paste0(gseaJar," GseaPreranked")
		}else{
			runCommand=paste0("java -Xmx8198m -cp ",gseaJar," xtools.gsea.GseaPreranked")
		}
		runCommand = paste0(runCommand, " -zip_report true -gmx ",gseaDb,"/",gseaCategory,	" -collapse false -nperm 1000 -rnk ",preRankedGeneFile," -scoring_scheme weighted -make_sets true -rpt_label '",gseaCategoryName,"' -plot_top_x 20 -set_max 500 -set_min 15 -out ", gesaResultDir)
		print(runCommand)
		system(runCommand)
	}

  resultDirSubs<-list.dirs(gesaResultDir,recursive=FALSE,full.names=TRUE)
  newResultDirSubs<-unlist(lapply(resultDirSubs, function(x) {
    newDir = gsub("\\.GseaPreranked.*", "", x)
    file.rename(x, newDir)
    return(newDir)
  }))
  
  dt<-data.frame(Folder=newResultDirSubs)
  dt$GseaCategory<-gsub(paste0(gesaResultDir,"/"), "", dt$Folder)
  write.csv(dt, file=paste0(basename(preRankedGeneFile),".gsea.csv"), row.names=F)
  
  if (makeReport) {
    library('rmarkdown')
    rmarkdown::render(gseaReportTemplt,output_file=paste0(basename(preRankedGeneFile),".gsea.html"),output_dir=resultDir)
  }

  return(dt)
}

preRankedGeneFileTable=parSampleFile1

if(!exists("makeReport")){
  makeReport=TRUE
}

preRankedGeneFileTable=read.delim(preRankedGeneFileTable,header=F,as.is=T)

alldt<-NULL
resultDir=getwd()

for (i in 1:nrow(preRankedGeneFileTable)) {
  preRankedGeneFile=preRankedGeneFileTable[i,1]
  
  if(!file.exists(preRankedGeneFile)){
    stop(paste0("File not exists: " + preRankedGeneFile))
  }
}
  
for (i in 1:nrow(preRankedGeneFileTable)) {
	preRankedGeneFile=preRankedGeneFileTable[i,1]
	
	compName=preRankedGeneFileTable[i,2]

  dt<-runGSEA(preRankedGeneFile,resultDir=resultDir,makeReport=makeReport,gseaJar=gseaJar,gseaDb=gseaDb,gseaCategories=gseaCategories)
  alldt<-rbind(alldt, dt)
}

write.csv(alldt, file="GSEA_result.csv", row.names=F)

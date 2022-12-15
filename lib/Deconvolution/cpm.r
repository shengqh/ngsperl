rm(list=ls())
outFile='discovery_cohort'
parSampleFile1=''
parSampleFile2=''
parSampleFile3=''
parFile1='/scratch/cqs/ravi_shah_projects/20221121_rnaseq_discovery_hg38/genetable/result/discovery_cohort.count'
parFile2='/data/cqs/references/deconvolution/tsp_v1_basisMatrix.txt'
parFile3=''


setwd('/scratch/cqs/ravi_shah_projects/20221121_rnaseq_discovery_hg38/deconvolution_cpm/result')

### Parameter setting end ###



library("edgeR")
library('data.table')


tsp<-fread(parFile2)

exp<-data.frame(fread(parFile1))
exp$Feature_gene_name<-gsub("-", ".", exp$Feature_gene_name)
feature_map<-unlist(split(exp$Feature_gene_name, exp$Feature))

# head(exp$Feature_gene_name)
# head(tsp$NAME)

# common<-exp[exp$Feature_gene_name %in% tsp$NAME,]
# table(common$Feature_gene_biotype)

counts<-exp[,c(8:ncol(exp))]

d <- DGEList(counts=counts)
d <- calcNormFactors(d, method="TMM")

tmmScale<-function(counts, libsize, tmm, mLplasma){
  scaleFac = tmm * libsize * mLplasma
  scaled <- sweep(counts, 2, scaleFac, FUN = '/')  * 10 ** 6
  return(scaled)
}

#since all plasma volumes are same, we will get identical result as edgeR cpm
ts<-tmmScale(counts, d$samples$lib.size, d$samples$norm.factors, 1)

#d_cpm <- cpm(d)
#head(d_cpm[,1:5])
#head(ts[,1:5])

ts$NAME<-feature_map[exp$Feature]
ts<-ts[ts$NAME %in% tsp$NAME,c("NAME", colnames(counts))]
write.csv(ts, paste0(outFile, ".cpm.csv"), row.names=F, quote=F)


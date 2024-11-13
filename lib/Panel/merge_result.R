rm(list=ls()) 
outFile='P11320'
parSampleFile1=''
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/h_cqs/yangj22/DavidHaas/P11320_MP_20240531/bwa_g4_refine_gatk4_SNV_12_panel_GenotypeGVCFs/result/P11320.g.vcf.gz'
parFile2='/nobackup/h_cqs/yangj22/DavidHaas/P11320_MP_20240531/bwa_g4_refine_gatk4_SNV_08_annovar/result/P11320/P11320.maf_filtered.annovar.final.tsv'
parFile3='/nobackup/h_cqs/yangj22/DavidHaas/P11320_MP_20240531/bwa_g4_refine_gatk4_SNV_03_GenotypeGVCFs/result/interval_for_gatk_GenotypeGVCFs_formatted.bed'
parFile4='/home/shengq2/program/collaborations/pipelines/panel_genes.txt'


setwd('/nobackup/h_cqs/yangj22/DavidHaas/P11320_MP_20240531/bwa_g4_refine_gatk4_SNV_13_panel_merge_result/result')

### Parameter setting end ###

###############################
####
#### 2024.09.27 Jing Yang
#### genotypes and phenotypes
###############################

library(data.table)

vcf_pre <- fread(parFile2, header=T, data.table=F)
print("VCF dim:")
dim(vcf_pre)

print("VCF FILTER:")
table(vcf_pre$FILTER)

first_sample_col=which(colnames(vcf_pre)=="FORMAT") + 1
vcf <- cbind(vcf_pre[,c("Chr","Start","avsnp150","Ref","Alt")],".",vcf_pre[,"FILTER"],vcf_pre$Gene.refGene,"GT",apply(vcf_pre[,first_sample_col:ncol(vcf_pre)],2,function(x) substr(x,1,12)))
colnames(vcf)[1:9] <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

gvcf_raw <- fread(parFile1, header=T, skip="CHROM", data.table=F)
print("GVCF dim:")
dim(gvcf_raw)

must_have_snps <- read.table(parFile3,header=F,sep="\t")
print("Must have SNPs dim:")
dim(must_have_snps)

gvcf_pre <- merge(must_have_snps,gvcf_raw,by.x=c("V1", "V3"),by.y=c("#CHROM", "POS"),all=T)
indels = which(gvcf_pre$V3 != gvcf_pre$V2 + 1)
print(paste("Number of indels:", length(indels)))

write.csv(gvcf_pre, file = "gvcf_pre_filling_indels.csv", quote = F, row.names = F)
indel=indels[1]
for(indel in indels){
  indel_start = gvcf_pre[indel,3]
  indel_end = gvcf_pre[indel,2]
  i = indel-1
  while(i > 0){
    if(gvcf_pre[i,1] < indel_start){
      break
    }
    if(is.na(gvcf_pre[i,3])){
      gvcf_pre[i,c(4:12)] <- gvcf_pre[indel,c(4:12)]
      gvcf_pre[i,3] <- gvcf_pre[i,2] - 1
    }
    i = i - 1
  }
}
write.csv(gvcf_pre, file = "gvcf_post_filling_indels.csv", quote = F, row.names = F)

print("GVCF dim:")
dim(gvcf_pre)

sample_column_index=which(colnames(gvcf_pre)=="FORMAT") + 1
gvcf <- cbind(gvcf_pre[,c("V1","V3","V4","REF","ALT")],".","PASS",gvcf_pre[,c("V11")],"GT",apply(gvcf_pre[,sample_column_index:ncol(gvcf_pre)],2,function(x) substr(x,1,12)))
colnames(gvcf)[1:9] <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

vcf_final <- unique(rbind(vcf,gvcf))
print("Final VCF dim:")
dim(vcf_final)

intgenes <- sort(readLines(parFile4))
all_files <- paste0(getwd(), "/", intgenes, "_annotation.txt")

# pharmcat <- read.table("/Users/jingyang/Dropbox/work/Other_projects/David_Haas_11320_20240530/pharmCat/pharmCat_output.txt",header = T,sep="\t")
pharmcat <- NULL
i=1
for (i in 1:length(intgenes)) {
  print(paste(i, "/", length(intgenes), ":", intgenes[i]))
  gene_output_file <- all_files[i]

  vcf_sub <- vcf_final[vcf_final$INFO %in% intgenes[i],]
  if (nrow(vcf_sub)!=0) {
    int_output <- matrix(NA,ncol(vcf_final)-9,nrow(vcf_sub))
    rownames(int_output) <- colnames(vcf_final)[10:ncol(vcf_final)]
    colnames(int_output) <- vcf_sub$ID
    for (k in 1:ncol(int_output)) {
      if (vcf_sub[k,"ID"] !="") {
        colnames(int_output)[k] <- vcf_sub[k,"ID"]
      } else {
        colnames(int_output)[k] <- paste0(vcf_sub[k,"#CHROM"],":",vcf_sub[k,"POS"])
      }
    }
    
    for (m in 1:nrow(vcf_sub)) {
      for(n in 10:ncol(vcf_sub)) {
        if (is.na(vcf_sub[m,rownames(int_output)[n-9]])) {
          output <- NA
        } else {
          tmp <- unlist(strsplit(vcf_sub[m,rownames(int_output)[n-9]],":"))
          output <- NA
          if (substr(tmp[1],1,1) == 0 & substr(tmp[1],3,3) == 0) {
            tmp2 <- unlist(strsplit(tmp[2],","))
            if (sum(as.numeric(tmp2))==0) {
              output <- "-/-"
            } else {
              output <- paste0(vcf_sub[m,"REF"],"/",vcf_sub[m,"REF"])
            }
          }
          if (substr(tmp[1],1,1) == 0 & substr(tmp[1],3,3) == 1 | substr(tmp[1],1,1) == 1 & substr(tmp[1],3,3) == 0) {
            output <- paste0(vcf_sub[m,"REF"],"/",vcf_sub[m,"ALT"])
          }  
          if (substr(tmp[1],1,1) == 0 & substr(tmp[1],3,3) == 2 | substr(tmp[1],1,1) == 2 & substr(tmp[1],3,3) == 0) {
            output <- paste0(vcf_sub[m,"REF"],"/",unlist(strsplit(vcf_sub[m,"ALT"],","))[2])
          }
          if (substr(tmp[1],1,1) == 1 & substr(tmp[1],3,3) == 1) {
            output <- paste0(vcf_sub[m,"ALT"],"/",vcf_sub[m,"ALT"])
          }
        }
        int_output[n-9,m] <- output
      } 
    }
    if (intgenes[i] %in% colnames(pharmcat)) {
      int_output <- data.frame(samples=rownames(int_output),pharmcat=pharmcat[,intgenes[i]],int_output)
      write.table(int_output,file = gene_output_file,quote=F,sep="\t",row.names = F)
    } else {
      int_output <- data.frame(samples=rownames(int_output),pharmcat=NA,int_output)
      write.table(int_output,file = gene_output_file,quote=F,sep="\t",row.names = F)
    }
  } else {
    int_output <- NULL
    write.table(int_output,file = gene_output_file,quote=F,sep="\t",row.names = F)
  }
}
writeLines(all_files, paste0(outFile, ".files.txt"))

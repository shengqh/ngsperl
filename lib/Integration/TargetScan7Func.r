library(data.table)

# Function 1 : get target mRNAs of candidate miRNAs
# miRNA_vec = vector of miRNAs
# targetscan_folder = folder with target scan database
# species = mmu or hsa
# prediction:
# Both_conserved: conserved sites of conserved miRNA families
# conserved: conserved and non conserved sites of conserved miRNA families
# non conserved: conserved and non conserved sites of non conserved miRNA families
# all: conserved and non conserved sites of conserved and non conserved miRNA families
# project_name: name that will be used to save txt

targetscan <- function(miRNA_vec, targetscan_folder, species=c("hsa","mmu"), prediction=c('Both_conserved', 'Conserved','Nonconserved','All'), save_targets= FALSE, project_name=NULL) {
  species<-match.arg(species)

  prediction<-match.arg(prediction)
  
  if (is.null(project_name)){
    project_name <- gsub("-", "", Sys.Date())
  }
  
  if(species == "hsa"){
    taxoid = "9606"
  }else if(species == "mmu"){
    taxoid = "10090"
  }
  
  family_file = paste0(targetscan_folder, "/", paste0(taxoid, "_miR_Family_Info.txt"))
  
  if(prediction == "Both_conserved"){
    data_file = paste0(targetscan_folder, "/", paste0(taxoid, "_Predicted_Targets_Info.default_predictions.txt"))
  }else{
    data_file = paste0(targetscan_folder, "/", paste0(taxoid, "_", prediction, "_Family_Info.txt"))
  }
  
  family <- fread(family_file, data.table = F, select = c(1,4))
  pred <- fread(data_file, data.table = F)
  
  miRNA_vec <- c(miRNA_vec)
  
  # get miRNA family 
  miRNA_family <- family[match(miRNA_vec,family$`MiRBase ID`), ]
  missing_miRNA<-setdiff(miRNA_vec,family$`MiRBase ID`)
  if (length(missing_miRNA) >0) {
    warning(paste0("
Not all miRNAs provided could match with the TargetScan family database, make sure you provide correct miRBase IDs. miRNAs not matched: ", paste0(missing_miRNA, collapse = ",")))
    fileConn<-file(paste0(project_name, "_miRNA_missing_in_database.txt"))
    writeLines(missing_miRNA, fileConn)
    close(fileConn)
  }
  
  # get miRNA targets
  miRNA_targets <- merge(miRNA_family, pred, by.x="miR family",by.y="miR Family", sort = FALSE)
  
  if(save_targets){
    write.table(miRNA_targets, paste0(project_name, "_miRNA_targets.txt"), sep="\t", quote = F)
  }
  
  return(miRNA_targets)
}

#test <- targetscan(as.vector(miRNA[,1]), targetscan_folder = "T:/Shared/Labs/Vickers Lab/MARS/TargetScan_7_2", species = "hsa", prediction = "Cons")


##############################################################################################################
# Function 2 : get target mRNAs of candidate miRNAs, match it with mRNA list. No fold change data needed.
  # miRNA = vector of miRNAs
  # mRNA = vector of mRNAs
  # species = mmu or hsa
  # prediction:
    # conserved: conserved and non conserved sites of conserved miRNA families
    # non conserved: conserved and non conserved sites of non conserved miRNA families
    # all: conserved and non conserved sites of conserved and non conserved miRNA families
  # save_targets: save txt file containin all targets of miRNAs (Irrespective of mRNA list)*
  # project_name: name that will be used to save txt
# *A txt file with the miRNA targtes that are also present in the provided mRNA vector is always saved

targetscan_mRNA <- function(miRNA_vec, mRNA_vec, targetscan_folder, species=c("hsa","mmu"), prediction=c('Conserved','Nonconserved','All'), save_targets= FALSE, project_name=NULL) {
  # get miRNA targets
  miRNA_targets <- targetscan(miRNA_vec, targetscan_folder, species, prediction, save_targets, project_name=project_name)
  
  #merge TargetScan targets with mRNA list
  miRNA_mRNA <- miRNA_targets[which(miRNA_targets$`Gene Symbol` %in% mRNA_vec),]
  
  if (save_targets) {
    write.table(miRNA_mRNA, paste(project_name, "_miRNA_mRNA_targets.txt", sep="_"), sep="\t", quote = F)
  }
  
  return(miRNA_mRNA)
}
# test <- targetscan_mRNA(miRNA, mRNA, species = "hsa", prediction = "conserved")


############################################################################################
# Function 3: Get target mRNAs of miRNAs, match up/down miRNAs with down/up mRNAs
  # miRNA = data.frame of DE miRNAs with FC
  # mRNA = data.frame of DE mRNAs with FC
  # species = mmu or hsa
  # prediction:
  #   conserved: conserved and non conserved sites of conserved miRNA families
  #   non conserved: conserved and non conserved sites of non conserved miRNA families
  #   all: conserved and non conserved sites of conserved and non conserved miRNA families
  # save_targets: save txt file containin all targets of miRNAs (Irrespective of mRNA list)*
  # project_name: name that will be used to save txt
  # *A txt file with the miRNA targtes that are also present in the provided mRNA vector is always saved

targetscan_DE <- function(miRNA, mRNA, targetscan_folder, species=c("hsa","mmu"), prediction=c('Conserved','Nonconserved','All'), save_targets= FALSE, project_name=NULL) {
  # get miRNA targets
  miRNA_targets <- targetscan(miRNA[,1], targetscan_folder, species, prediction, FALSE, project_name=project_name)
  
  miRNA_targets_foldchange <- merge(miRNA_targets, miRNA, by.x="MiRBase ID", by.y= "Feature_miRNA_name")
  
  #merge TargetScan targets with mRNA data.frame
  miRNA_mRNA_targets <- merge(miRNA_targets_foldchange, mRNA, by.x="Gene Symbol", by.y= "Feature_gene_name")
  
  if ( save_targets ) {
    write.table(miRNA_mRNA_targets, paste0(project_name, "_miRNA_mRNA_targets.txt"), sep="\t", quote = F)
  }
  
  #extract up/down miRNA and down/up mRNA. Save in separate files
  upmiRNA_dnmRNA <- miRNA_mRNA_targets[which(miRNA_mRNA_targets$log2FoldChange_miRNA > 0 & miRNA_mRNA_targets$log2FoldChange_mRNA < 0), ]
  upmiRNA_dnmRNA <<- upmiRNA_dnmRNA[,c("miR family", "MiRBase ID", "log2FoldChange_miRNA", "Gene Symbol", "Gene ID", "Transcript ID",
                                      "log2FoldChange_mRNA", "Species ID", "UTR start", "UTR end", "MSA start", "MSA end", "Seed match", "PCT")]
  
  write.table(upmiRNA_dnmRNA, paste0(project_name, "_TargetScan_UPmiRNA_DNmRNA.txt"), sep="\t", row.names = F, quote = F)

  dnmiRNA_upmRNA <- miRNA_mRNA_targets[which(miRNA_mRNA_targets$log2FoldChange_miRNA < 0 & miRNA_mRNA_targets$log2FoldChange_mRNA > 0), ]
  dnmiRNA_upmRNA <<- dnmiRNA_upmRNA[,c("miR family", "MiRBase ID", "log2FoldChange_miRNA", "Gene Symbol", "Gene ID", "Transcript ID",
                                      "log2FoldChange_mRNA", "Species ID", "UTR start", "UTR end", "MSA start", "MSA end", "Seed match", "PCT")]
  write.table(dnmiRNA_upmRNA, paste0(project_name, "_TargetScan_DNmiRNA_UPmRNA.txt"), sep="\t", row.names = F, quote = F)
}

#test <- targetscan_DE(miRNA, mRNA, species = "hsa", prediction = "all")



# MOFA functions

load_install <- function(package, method = "Default", install_package=package){
  if(!require(as.character(package),character.only = T)){
    
    if(method=="BiocM"){
      BiocManager::install(pkgs = install_package)
    }else{
      if(method=="Default"){
        install.packages(install_package)
      }else{
        stop("No valid method ('Default' for CRAM, 'BiocM' for BiocManager")
      }
    }
  }
  library(package,character.only = T)
}

create_mofa_object_long = function(named_list_of_dfs, feature_name_cols, sample_name_cols, value_name_cols){
  
  require("MOFA2")
  
  tog = as.data.frame(matrix(nrow = 0, ncol = 4))
  names(tog) = c("sample", "feature", "view", "value")
  for(i in 1:length(named_list_of_dfs)){
    the_obj = as.data.frame(named_list_of_dfs[[i]])
    cat(paste0("Head of features for ", names(named_list_of_dfs)[i], ": ", paste0(unlist(head(the_obj[,feature_name_cols[[i]]])), collapse = ", ")))
    cat(paste0("\nHead of samples for ", names(named_list_of_dfs)[i], ": ", paste0(unlist(head(the_obj[,sample_name_cols[[i]]])), collapse = ", ")))
    cat(paste0("\nHead of values for ", names(named_list_of_dfs)[i], ": ", paste0(unlist(head(the_obj[,value_name_cols[[i]]])), collapse = ", ")))
    cat("\n\n")
    
    to_add = as.data.frame(matrix(nrow = nrow(the_obj), ncol = 4))
    names(to_add) = c("sample", "feature", "view", "value")
    to_add$sample = the_obj[,sample_name_cols[[i]]]
    to_add$feature = the_obj[,feature_name_cols[[i]]]
    to_add$value = as.numeric(the_obj[,value_name_cols[[i]]])
    to_add$view = names(named_list_of_dfs)[i]
    
    tog = rbind(tog, to_add)
  }
  
  mofa_obj = create_mofa_from_df(tog)
  
  return(mofa_obj)
  
}

Uniprot_to_geneName = function(IDs, idmapping_file){
  uni_IDs = read.delim(idmapping_file, sep = "\t")
  uniProt_to_geneName = uni_IDs[which(uni_IDs$UniProtKB.ID=="Gene_Name"),]
  
  the_geneNames = uniProt_to_geneName$X1433B_HUMAN[match(IDs, uniProt_to_geneName$P31946)]
  
  return(the_geneNames)
  
}

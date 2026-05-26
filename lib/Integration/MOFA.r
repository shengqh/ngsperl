rm(list=ls()) 
outFile='PIMA'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20260525_PIMA_multiomics_reproduce/MOFA/result')

### Parameter setting end ###

source('MOFA_functions.R')

library(data.table)
library(dplyr)

myoptions_df=fread(parSampleFile1, header=F)

data_file_df=fread(parSampleFile2, header=F)
data_file_map=split(data_file_df$V1, data_file_df$V2)

def_df=fread(parSampleFile3, header=F)
mod_names=unique(def_df$V3)

get_cat_value = function(def_df, mod_name, cat_name){
  def_df |> 
    dplyr::filter(V3 == mod_name & V2 == cat_name) |> 
    dplyr::pull(V1)
}

mod1_name=mod_names[1]
mod1_file=data_file_map[[mod1_name]]
mod1_object_name=get_cat_value(def_df, mod1_name, "ObjectName")
mod1_sample_id_col=get_cat_value(def_df, mod1_name, "SampleColumn")
mod1_feature_name_col=get_cat_value(def_df, mod1_name, "FeatureColumn")
mod1_value_col=get_cat_value(def_df, mod1_name, "ValueColumn")

mod2_name=mod_names[2]
mod2_file=data_file_map[[mod2_name]]
mod2_object_name=get_cat_value(def_df, mod2_name, "ObjectName")
mod2_sample_id_col=get_cat_value(def_df, mod2_name, "SampleColumn")
mod2_feature_name_col=get_cat_value(def_df, mod2_name, "FeatureColumn")
mod2_value_col=get_cat_value(def_df, mod2_name, "ValueColumn")

load_rdata_env = function(path) {
  env = new.env(parent = baseenv())
  load(path, envir = env)
  env
}

require_object = function(mod_file, object_name, sample_id_col, feature_name_col, value_col) {
  env = load_rdata_env(mod_file)

  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Expected object '", object_name, "' in ", mod_file, call. = FALSE)
  }

  df = as.data.frame(get(object_name, envir = env, inherits = FALSE))

  if(any(!c(sample_id_col, feature_name_col, value_col) %in% colnames(df))){
    missed_cols = setdiff(c(sample_id_col, feature_name_col, value_col), colnames(df))
    stop(
      paste("Missed columns '", paste(missed_cols, collapse = ", "), "' in object '", object_name, "' in ", mod_file, " with colnames: ", paste(colnames(df), collapse = ", "), sep = ""),
      call. = FALSE
    )
  }

  df = df[, c(sample_id_col, feature_name_col, value_col)] |> 
    dplyr::rename(
      SampleID = all_of(sample_id_col),
      FeatureName = all_of(feature_name_col),
      Value = all_of(value_col)
    )
  return(df)
}

mod1_obj = as.data.frame(require_object(mod1_file, mod1_object_name, mod1_sample_id_col, mod1_feature_name_col, mod1_value_col))
  
mod2_obj = as.data.frame(require_object(mod2_file, mod2_object_name, mod2_sample_id_col, mod2_feature_name_col, mod2_value_col))

named_list_of_dfs = list(mod1_obj, mod2_obj)
names(named_list_of_dfs) = c(mod1_name, mod2_name)

message("Creating MOFA object")
mofa_obj = create_mofa_object_long(
  named_list_of_dfs = list(
    "Metabolites" = mod1_obj,
    "Proteomics" = mod2_obj
  ),
  feature_name_cols = list("FeatureName", "FeatureName"),
  sample_name_cols = list("SampleID", "SampleID"),
  value_name_cols = list("Value", "Value")
)

data_opts = get_default_data_options(mofa_obj)
model_opts = get_default_model_options(mofa_obj)
train_opts = get_default_training_options(mofa_obj)
train_opts$stochastic = FALSE

mofa_obj = prepare_mofa(
  object = mofa_obj,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

message("Running MOFA training")
output_hdf5=paste0(outFile, ".MOFA_model.hdf5")
trained_model = run_mofa(mofa_obj, output_hdf5, use_basilisk = FALSE)

message("MOFA model written to: ", output_hdf5)
invisible(trained_model)

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')

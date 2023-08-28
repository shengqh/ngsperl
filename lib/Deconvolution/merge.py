# merge together all the csv files, one for support vectors and the other for coefficients
# this script goes into each of the deconvolved sample directories
# and merges the coefficient files and the support vector files
# the best coefficient will then be chosen locally 
# Modified by Quanhu Sheng, based on https://github.com/sevahn/deconvolution/blob/master/deconvolve_cfrna_tutorial/merge_2.py

import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import logging

import os

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

def readFileMap(fileName):
  check_file_exists(fileName)

  result = OrderedDict()
  with open(fileName) as fh:
    for line in fh:
      filepath, name = line.strip().split('\t', 1)
      result[name] = filepath.strip()
  return(result)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="merge all the deconvolution results",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input file list', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file prefix", required=NotDEBUG)
args = parser.parse_args()

if DEBUG:
  listFile = "/nobackup/shah_lab/shengq2/20230628_rnaseq_discovery_hg38_validation/deconvolution_3_merge/result/discovery_cohort__fileList1.list"
  datasetName = "/nobackup/shah_lab/shengq2/20230628_rnaseq_discovery_hg38_validation/deconvolution_3_merge/result/discovery"
else:
  listFile = args.input
  datasetName = args.output

logger = logging.getLogger('merge')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

# read in samples
samps = pd.read_csv(listFile, sep = "\t", index_col = None, header = None) 
samps.columns = ["file_path", "sample_name"]
#samps.head()

sampNames = sorted(list(set(samps.sample_name.tolist())))
#sampNames

allSuppVec = pd.DataFrame()
allCoefs = pd.DataFrame()

# concatenate all the samples over all hyperparam combinations
for bioRep in sampNames:
  logger.info("Processing %s", bioRep)
  sub_df = samps[samps["sample_name"] == bioRep]
  sv_file_path = sub_df["file_path"][sub_df["file_path"].str.contains("supportVectors")].tolist()[0]
  coef_file_path = sub_df["file_path"][sub_df["file_path"].str.contains("Coefs")].tolist()[0]
  sampSV = pd.read_csv(sv_file_path, sep = ",", index_col = 0)
  sampCoef = pd.read_csv(coef_file_path, sep = ",", index_col = 0)
  allCoefs = pd.concat([sampCoef, allCoefs], ignore_index = False)
  allSuppVec = pd.concat([sampSV, allSuppVec], axis = "columns", ignore_index = False)

# pick the best hyperparam combination per sample
bestCoef = pd.DataFrame()

for samp in sampNames:
  allHyperThisSamp = [i for i in allCoefs.index if samp in i]
  hyperCombosThisSamp = allCoefs.loc[allHyperThisSamp]
  best = hyperCombosThisSamp[hyperCombosThisSamp['rmse'] == np.min(hyperCombosThisSamp['rmse'])]
  bestCoef = pd.concat([bestCoef, best])

performance = bestCoef.iloc[:,-2:].T
bestCoef = bestCoef.iloc[:,:-2]

# send negative coefs to zero
bestCoef[bestCoef < 0] = 0 # this is a samp x cells matrix
bestCoef = bestCoef.T

# normalize by total to get the fractions of cell type specific RNA
fracs = bestCoef.div(bestCoef.sum(axis = 0), axis = 1) 

# get the support vectors corresponding to the best coefficient pair
suppvecs = allSuppVec[bestCoef.columns.tolist()]

# strip the hyperparameter information so it's just the sample names
bestCoef.columns = [i.split("-NUSVR")[0] for i in bestCoef.columns.tolist()]
suppvecs.columns = [i.split("-NUSVR")[0] for i in suppvecs.columns.tolist()]
fracs.columns = [i.split("-NUSVR")[0] for i in fracs.columns.tolist()]

# write out the support vectors and the fractions
performance = performance.transpose()
performance["sample"] = bestCoef.columns
performance.to_csv(datasetName + ".performance.csv", sep = ",", header = True, index = True)

bestCoef.to_csv(datasetName + ".bestCoef.csv", sep = ",", header = True, index = True)
suppvecs.to_csv(datasetName + ".support_vectors.csv", sep = ",", header = True, index = False)
fracs.to_csv(datasetName + ".fractions.csv", sep = ",", header = True, index = True)

name_map = {
  "enterocyte of epithelium of large intestine/enterocyte of epithelium of small intestine/intestinal crypt stem cell of large intestine/large intestine goblet cell/mature enterocyte/paneth cell of epithelium of large intestine/small intestine goblet cell": "enterocyte/goblet/paneth of intestine",
  "club cell of prostate epithelium/hillock cell of prostate epithelium/hillock-club cell of prostate epithelium": "club/hillock of prostate epithelium",
  "intestinal enteroendocrine cell/paneth cell of epithelium of small intestine/transit amplifying cell of small intestine": "enteroendocrine/paneth/transit of intestine"
}  

def draw_figure(df, name_map, file_prefix, top_num = 5):
  top = list()
  for col in df.columns:
    cur_top = df.nlargest(top_num, col).index.tolist()
    top.extend(cur_top)

  top = sorted(list(set(top)))
  top_df = df.loc[top]
  top_renamed = top_df.rename(index = name_map)

  top_renamed.to_csv(file_prefix + ".top.csv", sep = ",", header = True, index = True)

  sns_plot = sns.heatmap(top_renamed, cmap = "RdBu_r", center = 0, xticklabels=True, yticklabels=True)
  sns_plot.figure.savefig(file_prefix + ".top.pdf", bbox_inches='tight')

draw_figure(fracs, name_map, datasetName + ".fractions")
draw_figure(bestCoef, name_map, datasetName + ".bestCoef")

logger.info("done.")
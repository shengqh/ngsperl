
import warnings
warnings.simplefilter(action='ignore')
import pycisTopic
pycisTopic.__version__
import scipy
import pandas as pd
import pickle 
from pycisTopic.cistopic_class import *
import argparse
import logging

DEBUG = False
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Create cisTopic object.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input fragment counts file', required=NOT_DEBUG)
parser.add_argument('-r', '--regions', action='store', nargs='?', help='Input regions file', required=NOT_DEBUG)
parser.add_argument('-c', '--cellnames', action='store', nargs='?', help='Input cell names file', required=NOT_DEBUG)
parser.add_argument('-m', '--meta', action='store', nargs='?', help='Input meta file', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output object pickle file", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="/nobackup/h_cqs/ramirema/maureen_gannon/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/scenic/10940_fragmentcounts.csv"
  args.regions="/nobackup/h_cqs/ramirema/maureen_gannon/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/scenic/10940_regions.txt"
  args.cellnames="/nobackup/h_cqs/ramirema/maureen_gannon/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/scenic/10940_cellnames.txt"
  args.meta="/nobackup/h_cqs/ramirema/maureen_gannon/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/scenic/10940_meta.csv"
  args.output="/nobackup/h_cqs/shengq2/temp/10940_snRNAseq_mmulatta_cisTopicObject.pkl"

logger = logging.getLogger('cisTopic')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("Read matrix from " + args.input)
count_matrix=scipy.io.mmread(args.input)
count_matrix=count_matrix.tocsr()

logger.info("Read regions from " + args.regions)
regions=pd.read_table(args.regions, header=None)
regions=regions[0].tolist()

logger.info("Read cell names from " + args.cellnames)
cellnames=pd.read_table(args.cellnames, header=None)
cellnames=cellnames[0].tolist()

logger.info("Create cistopic object")
cistopic_obj = create_cistopic_object(fragment_matrix=count_matrix, 
                                      cell_names=cellnames, 
                                      region_names=regions)
print(cistopic_obj.cell_data.head())

logger.info("Add cell data from " + args.meta)
cell_data = pd.read_csv(args.meta, index_col=0)
cell_data.rename(index = lambda x: x + "___cisTopic", inplace = True)
cistopic_obj.add_cell_data(cell_data)

print(cistopic_obj.cell_data.head())

logger.info("Save to " + args.output)
with open(args.output, 'wb') as f:
  pickle.dump(cistopic_obj, f)
  
logger.info("Done.")


import warnings
warnings.simplefilter(action='ignore')
import pycisTopic
pycisTopic.__version__
import pandas as pd
import pickle 
from pycisTopic.cistopic_class import *
import argparse
import logging
from pycisTopic.lda_models import run_cgs_models_mallet
import os

os.environ['MALLET_MEMORY'] = '200G'

DEBUG = False
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Run lda_models",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input object pickle file', required=NOT_DEBUG)
parser.add_argument('-m', '--mallet_path', action='store', nargs='?', help='Input mallet path', required=NOT_DEBUG)
parser.add_argument('-t', '--thread', action='store', nargs='?', type=int, help="Input threads", default=12)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output object pickle file", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/scenic_plus/T01_scenic_object/result/P10940.cisTopicObject.pkl"
  args.mallet_path="/data/cqs/ramirema/software/Mallet/bin/mallet"
  args.output="/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/scenic_plus/T01_scenic_object/result/P10940.lda_models.pkl"

logger = logging.getLogger('lda_models')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("Read pickle object from " + args.input)
with open(args.input, 'rb') as f:
  cistopic_obj = pickle.load(f)
  
#create temp folder "./tmp"
if not os.path.exists("./tmp"):
  os.makedirs("./tmp")

logger.info("Run models")
models=run_cgs_models_mallet(cistopic_obj, 
                             n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], 
                             n_cpu=args.thread, 
                             n_iter=500, 
                             random_state=555, 
                             alpha=50, 
                             alpha_by_topic=True, 
                             eta=0.1, 
                             eta_by_topic=False,
                             tmp_path="./tmp",
                             mallet_path=args.mallet_path)

logger.info("Save to " + args.output)
with open(args.output, 'wb') as f:
  pickle.dump(models, f)
  
logger.info("Done.")

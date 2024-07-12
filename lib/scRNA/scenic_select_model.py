
import warnings
warnings.simplefilter(action='ignore')
import pycisTopic
pycisTopic.__version__
import pickle 
from pycisTopic.cistopic_class import *
import argparse
import logging
from pycisTopic.lda_models import evaluate_models

DEBUG = False
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Run select model",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input object pickle file', required=NOT_DEBUG)
parser.add_argument('-m', '--models', action='store', nargs='?', help='Input models pickle file', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output object pickle file", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/scenic_plus/T01_scenic_object/result/P10940.cisTopicObject.pkl"
  args.models="/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/scenic_plus/T02_scenic_lda_models/result/P10940.lda_models.pkl"
  args.output="/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/scenic_plus/P10940.select_model.pkl"

logger = logging.getLogger('select_model')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("Read models from " + args.models)
with open(args.models, 'rb') as f:
  models = pickle.load(f)
  
logger.info("Select model")
model=evaluate_models(models,
                     select_model=None, 
                     return_model=True, 
                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics=False,
                     save= './model_selection.pdf')

logger.info("Read object from " + args.input)
with open(args.input, 'rb') as f:
  cistopic_obj = pickle.load(f)

logger.info("Add model to cisTopicObject")
cistopic_obj.add_LDA_model(model)

logger.info("Save to " + args.output)
with open(args.output, 'wb') as f:
  pickle.dump(cistopic_obj, f)
  
logger.info("Done.")

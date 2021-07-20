import os
import shutil
from zipfile import ZipFile

def extract_file(zipfile, extract_file, to_file, logger):
  logger.info("processing " + zipfile + " ...")
  with ZipFile(zipfile, 'r') as zipObj:
    listOfFiles = zipObj.namelist()
    f = [lf for lf in listOfFiles if lf.endswith(extract_file)][0]
    to_dir = os.path.dirname(to_file)
    ex_file = zipObj.extract(f, to_dir)
    shutil.copy(ex_file, to_file)

    ex_dir = os.path.dirname(ex_file)
    while ex_dir != to_dir:
      shutil.rmtree(ex_dir)
      ex_dir = os.path.dirname(ex_dir)
  logger.info("done")

if __name__ == "__main__":
    extract_file("/home/shengq2/program/ngsperl/data/shannon_vector.qzv","/data/metadata.tsv","/scratch/temp/shannon_vector.tsv")


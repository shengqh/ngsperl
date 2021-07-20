import os
from zipfile import ZipFile

def extract_file(zipfile, extract_file, to_file):
  zipfile = "/home/shengq2/program/ngsperl/data/shannon_vector.qzv"
  with ZipFile(zipfile, 'r') as zipObj:
    listOfFiles = zipObj.namelist()
    print(listOfFiles)

if __name__ == "__main__":
    extract_file("","","")



import sys
import os.path
import re

def get_name(filename):
  if(filename.endswith("FPKM")):
    filename = filename.replace(" FPKM", "_FPKM")
  
  if("Tag Count" in filename):
    filename = re.sub(" Tag Count.*", "_TagCount", filename)

  result = os.path.basename(filename)
  if result =="":
    result = filename

  return(result)

for matrix_file in sys.argv[1:]:
  print(f"Processing {matrix_file} ...")
  with open(matrix_file) as fin:
    data = fin.readlines()

  with open(matrix_file, 'w') as fout:
    header = data[0].rstrip().split('\t')
    snames = [get_name(x) if x.endswith('FPKM') or "Tag Count" in x else x for x in header[1:] ]
    new_header = "PeakID\t" + '\t'.join(snames) + "\n"
    fout.write(new_header)

    for line in data[1:]:
      fout.write(line)

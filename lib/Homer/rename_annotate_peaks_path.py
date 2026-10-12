
import sys
import os.path
import re

def get_name(filename):
  result = os.path.basename(filename)
  if result =="":
    result = filename

  result = result.replace('.tagDir', '').replace(" FPKM", "").replace("_FPKM", "")

  if (".NFR_filtered" in result):
    result = re.sub(r"\.NFR_filtered.*", "", result)

  if "Tag Count" in result or "_TagCount" in result:
    result = re.sub(r" Tag Count.*", "_TagCount", result)
    result = re.sub(r"_TagCount", "", result)

  return(result)

for matrix_file in sys.argv[1:]:
  print(f"Processing {matrix_file} ...")
  with open(matrix_file) as fin:
    data = fin.readlines()

  with open(matrix_file, 'w') as fout:
    header = data[0].rstrip().split('\t')
    print(header)
    snames = [get_name(x) if x.endswith('FPKM') or "Tag Count" in x or ".tagDir" in x or ".NFR_filtered" in x else x for x in header[1:] ]
    new_header = "PeakID\t" + '\t'.join(snames) + "\n"
    fout.write(new_header)

    for line in data[1:]:
      fout.write(line)


import sys
import os.path

def get_name(filenames_str):
  filenames = filenames_str.split('|')
  filenames = [os.path.basename(x).replace('.peaks.txt', '') for x in filenames]
  result = '|'.join(filenames)
  return(result)

peak_file = sys.argv[1]

with open(peak_file) as fin:
  data = fin.readlines()

with open(peak_file, 'w') as fout:
  for line in data:
    if line.startswith('#'):
      continue
    else:
      parts = line.split('\t')
      parts[6] = get_name(parts[6])
      fout.write('\t'.join(parts))

venn_file = sys.argv[2]
with open(venn_file) as fin:
  data = fin.readlines()

with open(venn_file, 'w') as fout:
  header = data[0].rstrip().split('\t')
  snames = [get_name(x) for x in header]
  new_header = '\t'.join(snames) + "\n"
  fout.write(new_header)

  for line in data[1:]:
    parts = line.rstrip().split('\t')
    parts[-1] = get_name(parts[-1])
    fout.write('\t'.join(parts) + '\n')

for matrix_file in sys.argv[3:]:
  with open(matrix_file) as fin:
    data = fin.readlines()

  with open(matrix_file, 'w') as fout:
    header = data[0].rstrip().split('\t')
    snames = [get_name(x) for x in header[1:]]
    new_header = "Sample\t" + '\t'.join(snames) + "\n"
    fout.write(new_header)

    for line in data[1:]:
      parts = line.rstrip().split('\t')
      parts[0] = get_name(parts[0])
      fout.write('\t'.join(parts) + '\n')

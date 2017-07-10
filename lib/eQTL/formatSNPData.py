import argparse

parser = argparse.ArgumentParser(description="format plink traw file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='plink traw format file', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', help='Output file', required=True)

args = parser.parse_args()

inputFile = args.input
outputFile = args.output

with open(inputFile, 'r') as f:
  with open(outputFile, 'w') as wr:
    headers = f.readline().strip().split('\t')
    nheaders = '\t'.join([h.split('_')[0] for h in headers])
    wr.write('%s\n' % nheaders)
    for line in f:
      wr.write(line)
     
    
    

import sys

inputfile = sys.argv[1]
outputfile = sys.argv[2]

with open(inputfile) as f:
  lis = [x.split() for x in f]

with open(outputfile, "w") as f:
  for x in zip(*lis):
    bFirst = True
    for y in x:
      if bFirst:
        f.write(y)
        bFirst = False
      else:
        f.write("\t" + y)
    f.write("\n")

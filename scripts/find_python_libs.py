import os
import os.path
from glob import glob

absfile=os.path.abspath(__file__)
PATH=os.path.dirname(os.path.dirname(absfile))
print(PATH)

result=[y for x in os.walk(PATH) for y in glob(os.path.join(x[0], '*.py'))]

ownlib=set()

#print(result)
libs=set()
for script in result:
  sfile=os.path.basename(script).replace(".py","")
  ownlib.add(sfile)
  with open(script, "rt") as fin:
    for line in fin:
      if "import" in line:
        sline=line.strip()
        if sline.startswith("import"):
          curlib=sline.replace("import","").strip()
          if " as " in sline:
            curlib = curlib.split(" as ")[0].strip()
          curlib=curlib.split(',')
          for l in curlib:
            libs.add(l.strip())
        elif sline.startswith("from"):
          curlib=sline.split('import')[0].replace("from","").split(',')
          for l in curlib:
            libs.add(l.strip())

alllibs=sorted([l for l in libs if l not in ownlib and not l.startswith("_")])
# alllibs=[l.replace("library(", "") for l in alllibs]
# alllibs=[l.replace(")", "") for l in alllibs]
# alllibs=[l.replace('"', "") for l in alllibs]
# alllibs=[l.replace("'", "") for l in alllibs]
# alllibs=[l for l in alllibs if not l.startswith('#')]
# alllibs=['BiocManager::install("' + l + '", update=FALSE, ask=FALSE, dependencies=TRUE)' for l in alllibs]

with open(os.path.dirname(absfile) + "/install_python_package.sh", "wt") as fout:
  for l in alllibs:
    fout.write(l + '\n')
    print(l)
